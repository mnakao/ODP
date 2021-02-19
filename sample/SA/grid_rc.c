// Relax constraints version
#include "common.h"

static int DISTANCE(const int u, const int v, const int height)
{
  int u_w = WIDTH(u,height);
  int u_h = HEIGHT(u,height);
  int v_w = WIDTH(v,height);
  int v_h = HEIGHT(v,height);
  return abs(u_w - v_w) + abs(u_h - v_h);
}

static int count_violating_edges(const int nodes, const int degree, const int adjacency[nodes][degree],
				 const int height, const int length)
{
  int c = 0;
  for(int i=0;i<nodes;i++)
    for(int j=0;j<degree;j++)
      if(DISTANCE(i, adjacency[i][j], height) > length)
	c++;
  
  return c;
}

static void print_help(char *argv)
{
  ERROR("%s -W width -H height -D degree -L length [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>] [-p <ratio of preprocessing>] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *width, int *height, int *degree, int *length, char *outfname, bool *enable_output,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, double *ratio, bool *enable_ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"W:H:D:L:o:s:n:w:c:p:A"))!=-1){
    switch(result){
    case 'W':
      *width = atoi(optarg);
      if(*width <= 0)
        ERROR("-W value > 0\n");
      break;
    case 'H':
      *height = atoi(optarg);
      if(*height <= 0)
        ERROR("-H value > 0\n");
      break;
    case 'D':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-D value > 0\n");
      break;
    case 'L':
      *length = atoi(optarg);
      if(*length <= 0)
        ERROR("-L value > 0\n");
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      strcpy(outfname, optarg);
      *enable_output = true;
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
    case 'p':
      *ratio =  atof(optarg);
      if(*ratio < 0 || *ratio > 1)
	ERROR("0.0 <= -p value <= 1.0\n");
        break;
    case 'A':
      *enable_ASPL_priority = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char outfname[MAX_FILENAME_LENGTH];
  bool enable_ASPL_priority = false, enable_output = false;
  int width = NOT_DEFINED, height = NOT_DEFINED, degree = NOT_DEFINED, length = NOT_DEFINED;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL, ratio = 0.5;

  set_args(argc, argv, &width, &height, &degree, &length, outfname, &enable_output,
	   &seed, &ncalcs, &max_temp, &min_temp, &ratio, &enable_ASPL_priority);
  int nodes = width * height;
  
  if(width == NOT_DEFINED || height == NOT_DEFINED || degree == NOT_DEFINED || length == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  
  printf("Width = %d, Height = %d, Degrees = %d, Length = %d\n",
	 width, height, degree, length);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  long pre_ncalcs = ncalcs * ratio;
  long post_ncalcs = ncalcs - pre_ncalcs;
  printf("Preprocessing ratio = %f (%ld)\n", ratio, pre_ncalcs);
  
  int lines = (nodes * degree)/2;
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  double create_time = get_time();
  ODP_Srand(seed);
  ODP_Generate_random_grid(width, height, degree, length+1, edge);
  create_time = get_time() - create_time;
  ODP_Conv_edge2adjacency_grid(width, height, lines, degree, edge, adjacency);
  
  ODP_Init_aspl_grid(width, height, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  // 1st
  ODP_Set_lbounds_grid(width, height, degree, length+1, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    ERROR("Unknown Error\n");
  }
  else{
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/pre_ncalcs);
    double temp = max_temp;
    int	interval = (pre_ncalcs < 100)? 1 : (pre_ncalcs)/100;
    printf("Ncalcs : Temp : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<pre_ncalcs;i++){
      if(i%interval == 0)
	printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      ODP_Mutate_adjacency_grid(width, height, degree, NULL, length+1, adjacency);
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
      
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, enable_ASPL_priority)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	ODP_Restore_adjacency_grid(adjacency);
      }
      temp *= cooling_rate;
    }
  }

  // 2nd
  ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/post_ncalcs);
    double temp = max_temp;
    int interval = (post_ncalcs < 100)? 1 : post_ncalcs/100;
    int tmp_length = length + 1, count = -1;
    bool once = true;
    printf("Ncalcs : Temp : Diameter : ASPL Gap\n");
    for(long i=0;i<post_ncalcs;i++){
      if(i%interval == 0)
        printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      ODP_Mutate_adjacency_grid(width, height, degree, NULL, length, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      if(count != 0){
	count = count_violating_edges(nodes, degree, adjacency, height, length);
	ASPL += count;
      }
      else{
	if(once){
	  best_diameter = diameter;
	  best_ASPL     = ASPL;
	  once          = false;
	}

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
      }

      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, enable_ASPL_priority)){
        current_diameter = diameter;
        current_ASPL     = ASPL;
      }
      else{
        ODP_Restore_adjacency_grid(adjacency);
      }
      temp *= cooling_rate;
    }
  }  
  
  sa_time = get_time() - sa_time;  
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_grid(width, height, degree, NULL, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f/%f sec. (Create Graph/SA)\n", create_time, sa_time);
  printf("ASPL Priority?  = %s\n", (enable_ASPL_priority)? "Yes" : "No");
  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");
  printf("Length          = %d\n", ODP_Get_length(lines, height, edge));

  if(enable_output){
    ODP_Write_edge_grid(lines, height, edge, outfname);
    printf("Generate ./%s\n", outfname);
  }
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
