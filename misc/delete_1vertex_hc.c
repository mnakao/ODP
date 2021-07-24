#include "common.h"

static void SWAP(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static void print_help(char *argv)
{
  ERROR("%s -f input.edges [-o output.edges] [-n ncalcs] [-s seed]\n", argv);
}

static void set_args(const int argc, char **argv, char **infname, char **outfname, long *ncalcs, int *seed)
{
  int result;
  while((result = getopt(argc,argv,"f:o:n:s:"))!=-1){
    switch(result){
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
    case 'n':
      *ncalcs = atol(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  int seed = 0, diameter = 0, new_diameter = 0, low_diameter = 0, new_low_diameter = 0, best_diameter = 0, current_diameter = 0;
  long ncalcs = 10000, sum = 0, new_sum = 0, best_sum = 0;
  double ASPL = 0, new_ASPL = 0, low_ASPL = 0, new_low_ASPL = 0, best_ASPL = 0, current_ASPL = 0;

  set_args(argc, argv, &infname, &outfname, &ncalcs, &seed);
  if(infname == NULL)
    ERROR("%s -f input.edges\n", argv[0]);

  ODP_Srand(seed); for(int i=0;i<100;i++) rand();
  int lines = ODP_Get_lines(infname);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  ODP_Read_edge_general(infname, edge);
  int nodes  = ODP_Get_nodes(lines, edge);
  int degree = ODP_Get_degree(nodes, lines, edge);
  if(degree%2 != 0) ERROR("degree must be an even number\n");
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);

  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();
  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Input file name = %s\n", infname);
  printf("Diameter        = %d\n", diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  
  int deleting_vertex =	get_random(nodes), new_nodes = nodes - 1, new_lines = lines - degree/2, partner[degree];
  int (*new_edge)[2]  = malloc(sizeof(int) * new_lines * 2);
  printf("---\n");
  printf("Deleting vertex %d\n", deleting_vertex);
  printf("---\n");

  for(int i=0,j=0,k=0;i<lines;i++){
    if(edge[i][0] != deleting_vertex && edge[i][1] != deleting_vertex){
      new_edge[j][0] = edge[i][0];
      new_edge[j][1] = edge[i][1];
      j++;
    }
    else if(edge[i][0] == deleting_vertex && edge[i][1] == deleting_vertex){
      ERROR("Loop is not supported\n");
    }
    else if(edge[i][0] == deleting_vertex){
      partner[k++] = edge[i][1];
    }
    else if(edge[i][1] == deleting_vertex){
      partner[k++] = edge[i][0];
    }
    else{
      ERROR("Someting Wrong 0\n");
    }
  }
    
  // Give randomness
  for(int i=0;i<degree*10;i++)
    SWAP(&partner[get_random(degree)], &partner[get_random(degree)]);

  memcpy(&new_edge[lines-degree][0], partner, sizeof(int)*degree);

  for(int i=0;i<new_lines;i++){
    if(new_edge[i][0] == deleting_vertex) ERROR("Someting Wrong 1\n");
    if(new_edge[i][0] >  deleting_vertex) new_edge[i][0]--;
    if(new_edge[i][1] >  deleting_vertex) new_edge[i][1]--;
  }

  for(int i=0;i<degree;i++){
    if(partner[i] == deleting_vertex)
      ERROR("Someting Wrong 2\n");
    else if(partner[i] > deleting_vertex)
      partner[i]--;
  }
  
  int (*new_adjacency)[degree] = malloc(sizeof(int) * new_nodes * degree);
  ODP_Conv_edge2adjacency_general(new_nodes, new_lines, degree, new_edge, new_adjacency);
  ODP_Set_lbounds_general(new_nodes, degree, &new_low_diameter, &new_low_ASPL);
  
  ODP_Init_aspl_general(new_nodes, degree, NULL);
  ODP_Set_aspl(new_adjacency, &new_diameter, &new_sum, &new_ASPL);
  best_diameter = current_diameter = new_diameter;
  best_sum      = new_sum;
  best_ASPL     = current_ASPL = new_ASPL;
  int (*best_adjacency)[degree] = malloc(sizeof(int) * new_nodes * degree);
  memcpy(best_adjacency, new_adjacency, sizeof(int) * new_nodes * degree);
  
  double hl_time = get_time();
  if(new_diameter == new_low_diameter && new_ASPL == new_low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    int tmp_partner[degree];
    long interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0)
        printf("%ld\t%d\t%f\n", i, best_diameter-new_low_diameter, best_ASPL-new_low_ASPL);

      memcpy(tmp_partner, partner, sizeof(int) * degree);
      SWAP(&partner[get_random(degree)], &partner[get_random(degree)]);
      memcpy(&new_edge[lines-degree][0], partner, sizeof(int)*degree);
      ODP_Conv_edge2adjacency_general(new_nodes, new_lines, degree, new_edge, new_adjacency);
      ODP_Set_aspl(new_adjacency, &new_diameter, &new_sum, &new_ASPL);

      if(new_diameter < best_diameter || (new_diameter == best_diameter && new_ASPL < best_ASPL)){
        best_diameter = new_diameter;
        best_sum      = new_sum;
        best_ASPL     = new_ASPL;
        memcpy(best_adjacency, new_adjacency, sizeof(int) * new_nodes * degree);
        if(new_diameter == new_low_diameter && new_ASPL == new_low_ASPL){
          printf("Find optimum solution\n");
          break;
        }
      }

      if(new_ASPL <= current_ASPL){
        current_diameter = new_diameter;
        current_ASPL     = new_ASPL;
      }
      else{
        memcpy(partner, tmp_partner, sizeof(int) * degree);
      }
    }
  }
  hl_time = get_time() - hl_time;
  ODP_Finalize_aspl();
  
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - new_low_diameter, best_diameter, new_low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)new_nodes*(new_nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - new_low_ASPL, best_ASPL, new_low_ASPL);
  printf("Time            = %f sec.\n", hl_time);

  if(outfname != NULL){
    printf("Output file name = %s\n", outfname);
    ODP_Conv_adjacency2edge_general(new_nodes, degree, NULL, best_adjacency, new_edge);
    ODP_Write_edge_general(new_lines, new_edge, outfname);
  }

  return 0;
}
