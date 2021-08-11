#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s -S symmetries [-N nodes] [-D degree] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-Y] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, int *symmetries, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *hill_climbing, bool *ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"N:D:S:f:o:s:n:w:c:YA"))!=-1){
    switch(result){
    case 'N':
      *nodes = atoi(optarg);
      if(*nodes <= 0)
        ERROR("-N value > 0\n");
      break;
    case 'D':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-D value > 0\n");
      break;
    case 'S':
      *symmetries = atoi(optarg);
      if(*symmetries <= 0)
        ERROR("-S value > 0\n");
      break;      
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      *outfname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*outfname, optarg);
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atol(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
      break;
    case 'w':
      *max_temp = atof(optarg);
      if(*max_temp <= 0)
        ERROR("-w value > 0\n");
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("MIN value > 0\n");
      break;
    case 'Y':
      *hill_climbing = true;
      break;      
    case 'A':
      *ASPL_priority = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static int NORM(int x, const int nodes)
{
  while(x < 0 || x >= nodes)
    x = (x < 0)? x + nodes : x - nodes;

  return x;
}

// return adjacency[v][d];
static int GLOBAL_ADJ(const int nodes, const int degree, const int symmetries,
                              const int (*adjacency)[degree], const int v, const int d)
{
  int based_nodes = nodes/symmetries;
  int n = adjacency[v%based_nodes][d] + (v/based_nodes)*based_nodes;
  return NORM(n, nodes);
}

static int LOCAL_INDEX(const int v, const int position, const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  return NORM(v - (position/based_nodes)*based_nodes, nodes);
}

static int get_index(const int u, const int v, const int u_d, const int nodes,
                     const int symmetries, const int degree, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ(nodes, degree, symmetries, adjacency, v, i) == u && i != u_d)
        return i;
  }
  else if(symmetries%2 == 0 && abs(u-v) == nodes/2){
    return u_d;
  }
  else{
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ(nodes, degree, symmetries, adjacency, v, i) == u)
        return i;
  }

  ERROR("Something Wrong ! [id=4]\n");
  return -1; // dummy
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool hill_climbing = false, ASPL_priority = false;
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, symmetries = 1, lines, (*edge)[2];
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ODP_Restore r;

  set_args(argc, argv, &nodes, &degree, &symmetries, &infname, &outfname, &seed,
           &ncalcs, &max_temp, &min_temp, &hill_climbing, &ASPL_priority);

  ODP_Srand(seed);
  if(infname){
    lines = ODP_Get_lines(infname);
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Read_edge_general(infname, edge);
    nodes  = ODP_Get_nodes(lines, edge);
    degree = ODP_Get_degree(nodes, lines, edge);
    if(nodes%symmetries != 0)
      ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);
  }
  else{
    if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
      print_help(argv[0]);
    else if(nodes%2 == 1 && degree%2 == 1)
      ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
    else if(nodes%symmetries != 0)
      ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);
    
    lines = (nodes * degree)/2;
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Generate_random_general_s(nodes, degree, symmetries, edge);
  }

  printf("Nodes = %d, Degrees = %d, Symmetries = %d\n", nodes, degree, symmetries);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  if(hill_climbing)
    printf("Method = Hill Climbing\n");
  else{
    printf("Method = Simulated Annealing\n");
    printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  }
  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);

  int based_nodes = nodes/symmetries;
  int (*adjacency)[degree] = malloc(sizeof(int)*based_nodes*degree);      // int adjacency[based_nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int)*based_nodes*degree); // int best_adjacency[based_nodes][degree];

  ODP_Conv_edge2adjacency_general_s(nodes, lines, degree, edge, symmetries, adjacency);

  // delete 1 edge
  int num_degrees[nodes], u, u_d, v, v_d;
  for(int i=0;i<nodes;i++) num_degrees[i] = degree;
  do{
    u   = get_random(nodes);
    u_d = get_random(degree);
    v   = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, u, u_d);
    v_d = get_index(u, v, u_d, nodes, symmetries, degree, adjacency);
  } while(u%based_nodes == v%based_nodes);

  int tmp = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, u, degree-1);
  adjacency[u%based_nodes][u_d] = LOCAL_INDEX(tmp, u, nodes, symmetries);
  adjacency[u%based_nodes][degree-1] = NOT_DEFINED;
  tmp = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, v, degree-1);
  adjacency[v%based_nodes][v_d] = LOCAL_INDEX(tmp, v, nodes, symmetries);
  adjacency[v%based_nodes][degree-1] = NOT_DEFINED;
  num_degrees[u%based_nodes]--;
  num_degrees[v%based_nodes]--;

  for(int i=based_nodes;i<nodes;i++)
    num_degrees[i] = num_degrees[i%based_nodes];

  /*
  adjacency[u][u_d]      = adjacency[u][degree-1];
  adjacency[u][degree-1] = NOT_DEFINED;
  adjacency[v][v_d]      = adjacency[v][degree-1];
  adjacency[v][degree-1] = NOT_DEFINED;
  num_degrees[u]--; num_degrees[v]--;
  */
  ODP_Init_aspl_general_s(nodes, degree, num_degrees, symmetries); 
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL); 

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*based_nodes*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = (hill_climbing)? 0 : pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    if(hill_climbing)
      printf("Ncalcs : Best ASPL Gap ( Dia. )\n");
    else
      printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0){
        if(hill_climbing)
          printf("%ld\t%f ( %d )\n", i,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
        else
          printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
                 current_ASPL-low_ASPL, current_diameter-low_diameter,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
      }

      ODP_Mutate_adjacency_general_s(nodes, degree, num_degrees, symmetries, &r, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*based_nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, hill_climbing, ASPL_priority, symmetries)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	ODP_Restore_adjacency_general(r, adjacency);
      }
      temp *= cooling_rate;
    }
  }

  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  for(int i=0;i<lines;i++)
    edge[i][0] = edge[i][1] = -1;
  ODP_Conv_adjacency2edge_general_s(nodes, degree, num_degrees, best_adjacency, symmetries, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  //  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  //  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");

  lines = nodes * degree / 2 - symmetries;
  if(outfname)
    ODP_Write_edge_general(lines, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
