#include "common.h"

static void SWAP(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static void print_help(char *argv)
{
  ERROR("%s -f input.edges [-o output.edges] [-s seed]\n", argv);
}

static void set_args(const int argc, char **argv, char **infname, char **outfname, int *seed)
{
  int result;
  while((result = getopt(argc,argv,"f:o:s:"))!=-1){
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
  int seed = 0, diameter = 0, new_diameter = 0, low_diameter = 0, new_low_diameter = 0;
  long sum = 0, new_sum = 0;
  double ASPL = 0, new_ASPL = 0, low_ASPL = 0, new_low_ASPL = 0;

  set_args(argc, argv, &infname, &outfname, &seed);
  if(infname == NULL)
    ERROR("%s -f input.edges\n", argv[0]);

  ODP_Srand(seed);
  int lines = ODP_Get_lines(infname);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  ODP_Read_edge_general(infname, edge);
  int nodes  = ODP_Get_nodes(lines, edge);
  int degree = ODP_Get_degree(nodes, lines, edge);
  if(degree%2 != 0) ERROR("degree must be an even number\n");
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Input file name = %s\n", infname);

  for(int i=0;i<100;i++) rand();
  int deleting_vertex =	get_random(nodes);
  int new_nodes       = nodes - 1;
  int new_lines       = lines - degree/2;
  int (*new_edge)[2]  = malloc(sizeof(int) * new_lines * 2);
  
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);
  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();
  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  
  printf("---\n");
  printf("Deleting vertex %d\n", deleting_vertex);
  printf("---\n");
  
  int partner[degree];
  int j = 0, k = 0;
  for(int i=0;i<lines;i++){
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
  
  if(j != lines-degree)
    ERROR("Someting Wrong 1\n");
  else if(k != degree)
    ERROR("Someting Wrong 2\n");
  
  // Give randomness
  for(int i=0;i<degree*10;i++)
    SWAP(&partner[get_random(degree)], &partner[get_random(degree)]);

  for(int i=lines-degree,j=0;i<new_lines;i++,j+=2){
    new_edge[i][0] = partner[j  ];
    new_edge[i][1] = partner[j+1];
  }

  for(int i=0;i<new_lines;i++){
    if(new_edge[i][0] == deleting_vertex) ERROR("Someting Wrong 3\n");
    if(new_edge[i][0] >  deleting_vertex) new_edge[i][0]--;
    if(new_edge[i][1] >  deleting_vertex) new_edge[i][1]--;
  }
  
  if(outfname == NULL){
    ODP_Print_edge_general(new_lines, new_edge);
  }
  else{
    printf("Output file name = %s\n", outfname);
    ODP_Write_edge_general(new_lines, new_edge, outfname);
  }

  int (*new_adjacency)[degree] = malloc(sizeof(int) * new_nodes * degree);
  ODP_Conv_edge2adjacency_general(new_nodes, new_lines, degree, new_edge, new_adjacency);
  ODP_Init_aspl_general(new_nodes, degree, NULL);
  ODP_Set_aspl(new_adjacency, &new_diameter, &new_sum, &new_ASPL);
  ODP_Finalize_aspl();
  ODP_Set_lbounds_general(new_nodes, degree, &new_low_diameter, &new_low_ASPL);
  printf("Diameter     = %d\n", new_diameter);
  printf("Diameter Gap = %d (%d - %d)\n", new_diameter - new_low_diameter, new_diameter, new_low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", new_ASPL, new_sum, (double)new_nodes*(new_nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", new_ASPL - new_low_ASPL, new_ASPL, new_low_ASPL);

  return 0;
}
