#define METRIC_MAX 16

typedef enum metric_type{
  metric_int,
  metric_double,
  metric_string
} metric_type;

static int metric_ticker = 0;
static metric_type types[METRIC_MAX];
static int ints[METRIC_MAX];
static int doubles[METRIC_MAX];
static char *strings[METRIC_MAX];
static char *keys[METRIC_MAX];

void metric_load_int(char *key, int value){
  if (metric_ticker >= METRIC_MAX) {
    printf("ReproBLAS error: too many metrics\n");
    exit(125);
  }
  keys[metric_ticker] = key;
  types[metric_ticker] = metric_int;
  ints[metric_ticker] = value;
  metric_ticker++;
}

void metric_load_double(char *key, double value){
  if (metric_ticker >= METRIC_MAX) {
    printf("ReproBLAS error: too many metrics\n");
    exit(125);
  }
  keys[metric_ticker] = key;
  types[metric_ticker] = metric_double;
  doubles[metric_ticker] = value;
  metric_ticker++;
}

void metric_load_string(char *key, char *value){
  if (metric_ticker >= METRIC_MAX) {
    printf("ReproBLAS error: too many metrics\n");
    exit(125);
  }
  keys[metric_ticker] = key;
  types[metric_ticker] = metric_string;
  strings[metric_ticker] = value;
  metric_ticker++;
}

void metric_dump(){
  int i;
  comma = ","
  printf("{\n");
  for(i = 0; i < metric_ticker; i++){
    if (i == metric_ticker - 1){
      comma = "";
    }
    switch(types[i]){
      case metric_int:
        printf("\t%s: %d%s\n", keys[i], ints[i], comma);
        break;
      case metric_double:
        printf("\t%s: %g%s\n", keys[i], doubles[i], comma);
        break;
      case metric_string:
        printf("\t%s: %g%s\n", keys[i], strings[i], comma);
        break;
    }
  }
  printf("}\n");
}
#endif
