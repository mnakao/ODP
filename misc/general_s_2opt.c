#include "common.h"
static int _u[2], _u_d[2], _v[2], _v_d[2], _nodes, _symmetries, _based_nodes, _degree;
static bool check_rotated_edges_overlap_general(const int u0, const int v0, const int u1, const int v1,
                                                const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  int diff0 = (u0 > v0)? v0 - u0 + nodes : v0 - u0;
  int diff1 = (u1 > v1)? v1 - u1 + nodes : v1 - u1;
  int diff2 = (u1 < v1)? u1 - v1 + nodes : u1 - v1;

  if(diff0 == diff1 && u0%based_nodes == u1%based_nodes)
    return true;
  else if(diff0 == diff2 && u0%based_nodes == v1%based_nodes)
    return true;
  else
    return false;
}

static void SWAP(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static bool IS_DIAMETER_GENERAL(const int u, const int v, const int nodes, const int symmetries)
{
  if(symmetries%2 != 0 || abs(u-v) != nodes/2)
    return false;
  else
    return true;
}

static int NORM(int x, const int nodes)
{
  while(x < 0 || x >= nodes)
    x = (x < 0)? x + nodes : x - nodes;

  return x;
}

static int LOCAL_INDEX_GENERAL(const int v, const int position, const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  return NORM(v - (position/based_nodes)*based_nodes, nodes);
}

static bool check_multiple_edges_general_s(const int u, const int u_d, const int v, const int nodes, const int degree,
                                           const int symmetries, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  for(int i=0;i<degree;i++)
    if(i!=u_d && adjacency[u][i] == v)
      return false;

  return true;
}

// return adjacency[v][d];
static int GLOBAL_ADJ_GENERAL(const int nodes, const int degree, const int symmetries,
                              const int (*adjacency)[degree], const int v, const int d)
{
  int based_nodes = nodes/symmetries;
  int n = adjacency[v%based_nodes][d] + (v/based_nodes)*based_nodes;
  return NORM(n, nodes);
}

static int get_degree_index_general(const int u, const int v, const int u_d, const int nodes,
                                    const int symmetries, const int degree, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u && i != u_d)
        return i;
  }
  else if(symmetries%2 == 0 && abs(u-v) == nodes/2){
    return u_d;
  }
  else{
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u)
        return i;
  }

  ERROR("Something Wrong ! [id=4]\n");
  return -1; // dummy
}

void restore_adjacency_general(int (*adjacency)[_degree])
{
  for(int i=0;i<2;i++){
    adjacency[_u[i]%_based_nodes][_u_d[i]] = LOCAL_INDEX_GENERAL(_v[i], _u[i], _nodes, _symmetries);
    adjacency[_v[i]%_based_nodes][_v_d[i]] = LOCAL_INDEX_GENERAL(_u[i], _v[i], _nodes, _symmetries);
  }
}

static void backup_restore_adjacency(const int u[2], const int u_d[2], const int v[2], const int v_d[2], const int nodes, const int symmetries,
                                     const int degree)
{
  _nodes = nodes;
  _symmetries = symmetries;
  _based_nodes = nodes/symmetries;
  _degree = degree;
  
  for(int i=0;i<2;i++){
    _u[i]   = u[i];
    _v[i]   = v[i];
    _u_d[i] = u_d[i];
    _v_d[i] = v_d[i];
  }
}

static bool mutate_adjacency_1opt_general_s(const int u, const int u_d, const int nodes, const int degree,
                                            const int symmetries, int adjacency[nodes/symmetries][degree])
{
  int based_nodes = nodes/symmetries;
  int v = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, (const int (*)[degree])adjacency, u, u_d);
  if(IS_DIAMETER_GENERAL(u, v, nodes, symmetries)) return false;
  int v_d = get_degree_index_general(u, v, u_d, nodes, symmetries, degree, (const int (*)[degree])adjacency);
  int new_v;
  if((u-v)%based_nodes == 0){
    if(symmetries <= 4) return false;
    while(1){
      new_v = v + based_nodes * get_random(symmetries);
      if(u == new_v || v == new_v) continue;
      else if(NORM(u-v, nodes) == NORM(new_v-u, nodes)) continue;
      else if(nodes%2 == 1) break;
      else if(/* nodes%2 == 0 && */ (u-new_v)%(nodes/2) != 0) break;
    }
  }
  else{
    int rnd   = (symmetries%2 == 1)? get_random(symmetries-1) : get_random(symmetries);
    new_v = (rnd != symmetries-1)? v + based_nodes*(rnd+1) : u + based_nodes*(symmetries/2);
    //  int rnd   = get_random(symmetries-1);
    //  new_v = v + based_nodes*(rnd+1);
  }
  int new_u = v - (new_v - u);
  int tmp[2] = {LOCAL_INDEX_GENERAL(new_v, u, nodes, symmetries), LOCAL_INDEX_GENERAL(new_u, v, nodes, symmetries)};
  if(!check_multiple_edges_general_s(u%based_nodes, u_d, tmp[0], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
     !check_multiple_edges_general_s(v%based_nodes, v_d, tmp[1], nodes, degree, symmetries, (const int (*)[degree])adjacency))
    return false;
  
  adjacency[u%based_nodes][u_d] = tmp[0];
  adjacency[v%based_nodes][v_d] = tmp[1];
  return true;
}

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
  ODP_Init_aspl_general_s(nodes, degree, NULL, symmetries); 
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

      while(1){
        int u[2], v[2], u_d[2], v_d[2];
        while(1){
          u[0] = get_random(nodes);
          u[1] = get_random(nodes);
          if(u[0] == u[1]) continue;

          u_d[0] = get_random(degree);
          v[0] = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, (const int (*)[degree])adjacency, u[0], u_d[0]);
          if(v[0] == u[1]) continue;
          
          u_d[1] = get_random(degree);
          v[1] = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, (const int (*)[degree])adjacency, u[1], u_d[1]);
          if(v[1] == u[0] || v[0] == v[1]) continue;
          break;
        }
        v_d[0] = get_degree_index_general(u[0], v[0], u_d[0], nodes, symmetries, degree, (const int (*)[degree])adjacency);
        v_d[1] = get_degree_index_general(u[1], v[1], u_d[1], nodes, symmetries, degree, (const int (*)[degree])adjacency);
        backup_restore_adjacency(u, u_d, v, v_d, nodes, symmetries, degree);

        if(IS_DIAMETER_GENERAL(u[0], v[0], nodes, symmetries) && IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
          if((u[0] - u[1])%based_nodes == 0)
            continue;

          int tmp[2];
          if(get_random(2)){
            tmp[0] = LOCAL_INDEX_GENERAL(u[1], u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(u[0], u[1], nodes, symmetries);
          }
          else{
            tmp[0] = LOCAL_INDEX_GENERAL(v[1], u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(v[0], u[1], nodes, symmetries);
          }
          if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
             !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, symmetries, (const int (*)[degree])adjacency))
            continue;

          adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
          adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
          break;
        }
        else if(IS_DIAMETER_GENERAL(u[0], v[0], nodes, symmetries) || IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
          if(IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
            SWAP(&u[0], &u[1]); SWAP(&u_d[0], &u_d[1]);
            SWAP(&v[0], &v[1]); SWAP(&v_d[0], &v_d[1]);
          }

          int opposite = nodes/2, rnd = get_random(4), tmp[4];
          if(rnd == 0){ // u[0]--v[1], u[1]--u[1]', v[0]--v[1]'
            tmp[0] = LOCAL_INDEX_GENERAL(v[1],          u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[1], nodes, symmetries);
            tmp[2] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[0], nodes, symmetries);
            tmp[3] = LOCAL_INDEX_GENERAL(u[0],          v[1], nodes, symmetries);
          }
          else if(rnd == 1){ // u[0]--v[1]', v[0]--v[1], u[1]--u[1]'
            tmp[0] = LOCAL_INDEX_GENERAL(v[1]+opposite, u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[1], nodes, symmetries);
            tmp[2] = LOCAL_INDEX_GENERAL(v[1],          v[0], nodes, symmetries);
            tmp[3] = LOCAL_INDEX_GENERAL(v[0],          v[1], nodes, symmetries);
          }
          else if(rnd == 2){ // u[0]--u[1], v[0]--u[1]', v[1]--v[1]'
            tmp[0] = LOCAL_INDEX_GENERAL(u[1],          u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(u[0],          u[1], nodes, symmetries);
            tmp[2] = LOCAL_INDEX_GENERAL(u[1]+opposite, v[0], nodes, symmetries);
            tmp[3] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[1], nodes, symmetries);
          }
          else if(rnd == 3){ // u[0]--u[1]', u[1]--v[0], v[1]--v[1]'
            tmp[0] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[0], nodes, symmetries);
            tmp[1] = LOCAL_INDEX_GENERAL(v[0],          u[1], nodes, symmetries);
            tmp[2] = LOCAL_INDEX_GENERAL(u[1],          v[0], nodes, symmetries);
            tmp[3] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[1], nodes, symmetries);
          }
          
          if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
             !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
             !check_multiple_edges_general_s(v[0]%based_nodes, v_d[0], tmp[2], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
             !check_multiple_edges_general_s(v[1]%based_nodes, v_d[1], tmp[3], nodes, degree, symmetries, (const int (*)[degree])adjacency))
            continue;
          
          adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
          adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
          adjacency[v[0]%based_nodes][v_d[0]] = tmp[2];
          adjacency[v[1]%based_nodes][v_d[1]] = tmp[3];
          break;
        }

        // Two selected edges are symmetrical
        if(check_rotated_edges_overlap_general(u[0], v[0], u[1], v[1], nodes, symmetries)){
          if(mutate_adjacency_1opt_general_s(u[0], u_d[0], nodes, degree, symmetries, adjacency)){
            break;
          }
          else{
            continue;
          }
        }
        
        int tmp[4];
        if(get_random(2)){ // u[0]--v[1], v[0]--u[1]
          if(IS_DIAMETER_GENERAL(u[0], v[1], nodes, symmetries) || IS_DIAMETER_GENERAL(v[0], u[1], nodes, symmetries))
            continue;
          else if(check_rotated_edges_overlap_general(u[0], v[1], v[0], u[1], nodes, symmetries))
            continue;
          tmp[0] = LOCAL_INDEX_GENERAL(v[1], u[0], nodes, symmetries);
          tmp[1] = LOCAL_INDEX_GENERAL(v[0], u[1], nodes, symmetries);
          tmp[2] = LOCAL_INDEX_GENERAL(u[1], v[0], nodes, symmetries);
          tmp[3] = LOCAL_INDEX_GENERAL(u[0], v[1], nodes, symmetries);
        }
        else{ // u[0]--u[1], v[0]--v[1]
          if(IS_DIAMETER_GENERAL(u[0], u[1], nodes, symmetries) || IS_DIAMETER_GENERAL(v[0], v[1], nodes, symmetries))
            continue;
          else if(check_rotated_edges_overlap_general(u[0], u[1], v[0], v[1], nodes, symmetries))
            continue;
          tmp[0] = LOCAL_INDEX_GENERAL(u[1], u[0], nodes, symmetries);
          tmp[1] = LOCAL_INDEX_GENERAL(u[0], u[1], nodes, symmetries);
          tmp[2] = LOCAL_INDEX_GENERAL(v[1], v[0], nodes, symmetries);
          tmp[3] = LOCAL_INDEX_GENERAL(v[0], v[1], nodes, symmetries);
        }

        if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
           !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
           !check_multiple_edges_general_s(v[0]%based_nodes, v_d[0], tmp[2], nodes, degree, symmetries, (const int (*)[degree])adjacency) ||
           !check_multiple_edges_general_s(v[1]%based_nodes, v_d[1], tmp[3], nodes, degree, symmetries, (const int (*)[degree])adjacency))
          continue;

        adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
        adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
        adjacency[v[0]%based_nodes][v_d[0]] = tmp[2];
        adjacency[v[1]%based_nodes][v_d[1]] = tmp[3];
        break;
      }

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
        restore_adjacency_general(adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_general_s(nodes, degree, NULL, best_adjacency, symmetries, edge);
  
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
