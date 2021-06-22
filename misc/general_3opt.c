#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s [-N nodes] [-D degree] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-Y] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *hill_climbing, bool *ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"N:D:f:o:s:n:w:c:YA"))!=-1){
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

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool hill_climbing = false, ASPL_priority = false;
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, lines, (*edge)[2];
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, &infname, &outfname, &seed,
	   &ncalcs, &max_temp, &min_temp, &hill_climbing, &ASPL_priority);
  
  ODP_Srand(seed);
  if(infname){
    lines = ODP_Get_lines(infname);
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Read_edge_general(infname, edge);
    nodes  = ODP_Get_nodes(lines, edge);
    degree = ODP_Get_degree(nodes, lines, edge);
  }
  else{
    if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
      print_help(argv[0]);
    else if(nodes%2 == 1 && degree%2 == 1)
      ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);

    lines = (nodes * degree)/2;
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Generate_random_general(nodes, degree, edge);
  }
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
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

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = (hill_climbing)? 0 : pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    long j = 0;
    if(hill_climbing)
      printf("Ncalcs : Best ASPL Gap ( Dia. )\n");
    else
      printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    for(long i=0;i<ncalcs;i++){
      if(i/interval == j){
        j++;
        if(hill_climbing)
          printf("%ld\t%f ( %d )\n", i, 
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
        else
          printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
                 current_ASPL-low_ASPL, current_diameter-low_diameter,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
      }

      int u[3], v[3], u_d[3], v_d[3], new_u[3], new_v[3];
      while(1){
        while(1){
          u[0] = get_random(nodes);
          u[1] = get_random(nodes);
          if(u[0] == u[1]) continue;
          u[2] = get_random(nodes);
          if(u[2] == u[0] || u[2] == u[1]) continue;

          u_d[0] = get_random(degree);
          v[0]   = adjacency[u[0]][u_d[0]];
          if(v[0] == u[1] || v[0] == u[2]) continue;
          
          u_d[1] = get_random(degree);
          v[1]   = adjacency[u[1]][u_d[1]];
          if(v[1] == u[0] || v[1] == v[0] || v[1] == u[2]) continue;

          u_d[2] = get_random(degree);
          v[2]   = adjacency[u[2]][u_d[2]];
          if(v[2] == u[0] || v[2] == u[1] || v[2] == v[0] || v[2] == v[1]) continue;
          break;
        }

        v_d[0] = get_degree_index(u[0], v[0], u_d[0], nodes, degree, (const int (*)[degree])adjacency);
        v_d[1] = get_degree_index(u[1], v[1], u_d[1], nodes, degree, (const int (*)[degree])adjacency);
        v_d[2] = get_degree_index(u[2], v[2], u_d[2], nodes, degree, (const int (*)[degree])adjacency);
        backup_adjacency(3, u, u_d, v, v_d);

        // 3-opt
        int r = get_random(8);
        if(r == 0){      // v[0]--v[1], u[1]--u[2], u[0]--v[2]
          new_u[0] = v[1]; new_u[1] = v[0]; new_u[2] = u[0];
          new_v[0] = v[2]; new_v[1] = u[2]; new_v[2] = u[1];
        }
        else if(r == 1){ // v[0]--v[1], u[1]--v[2], u[0]--u[2]
          new_u[0] = v[1]; new_u[1] = v[0]; new_u[2] = u[1];
      	  new_v[0] = u[2]; new_v[1] = v[2]; new_v[2] = u[0];
        }
        else if(r == 2){ // v[0]--u[1], v[1]--v[2], u[2]--u[0]
          new_u[0] = u[1]; new_u[1] = v[2]; new_u[2] = v[1];
          new_v[0] = u[2]; new_v[1] = v[0]; new_v[2] = u[0];
        }
        else if(r == 3){ // v[0]--u[1], v[1]--u[2], v[2]--u[0]
          new_u[0] = u[1]; new_u[1] = u[2]; new_u[2] = u[0];
          new_v[0] = v[2]; new_v[1] = v[0]; new_v[2] = v[1];
        }
        else if(r == 4){ // v[0]--v[2], u[0]--v[1], u[1]--u[2]
          new_u[0] = v[2]; new_u[1] = u[0]; new_u[2] = v[0];
          new_v[0] = v[1]; new_v[1] = u[2]; new_v[2] = u[1];
        }
        else if(r == 5){ // v[0]--v[2], u[0]--u[1], v[1]--u[2]
          new_u[0] = v[2]; new_u[1] = u[2]; new_u[2] = v[0];
          new_v[0] = u[1]; new_v[1] = u[0]; new_v[2] = v[1];
        }
        else if(r == 6){ // v[0]--u[2], u[0]--v[1], u[1]--v[2]
          new_u[0] = u[2]; new_u[1] = u[0]; new_u[2] = u[1];
          new_v[0] = v[1]; new_v[1] = v[2]; new_v[2] = v[0];
        }
        else if(r == 7){ // v[0]--u[2], u[0]--u[1], v[1]--v[2]
          new_u[0] = u[2]; new_u[1] = v[2]; new_u[2] = v[1];
          new_v[0] = u[1]; new_v[1] = u[0]; new_v[2] = v[0];
        }

        if(!check_multiple_edges(u[0], u_d[0], new_v[0], degree, (const int (*)[degree])adjacency) ||
           !check_multiple_edges(u[1], u_d[1], new_v[1], degree, (const int (*)[degree])adjacency) ||
           !check_multiple_edges(u[2], u_d[2], new_v[2], degree, (const int (*)[degree])adjacency) ||
           !check_multiple_edges(v[0], v_d[0], new_u[0], degree, (const int (*)[degree])adjacency) ||
           !check_multiple_edges(v[1], v_d[1], new_u[1], degree, (const int (*)[degree])adjacency) ||
           !check_multiple_edges(v[2], v_d[2], new_u[2], degree, (const int (*)[degree])adjacency))
          continue;

        adjacency[u[0]][u_d[0]] = new_v[0];
        adjacency[u[1]][u_d[1]] = new_v[1];
        adjacency[u[2]][u_d[2]] = new_v[2];
        adjacency[v[0]][v_d[0]] = new_u[0];
        adjacency[v[1]][v_d[1]] = new_u[1];
        adjacency[v[2]][v_d[2]] = new_u[2];
        break;       
      }
        
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, hill_climbing, ASPL_priority)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
        undo_adjacency(3, degree, adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_general(nodes, degree, NULL, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  //  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  //  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");

  if(outfname)
    ODP_Write_edge_general(lines, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
