#define METRIC_MAX 16

typedef enum metric_type{
  metric_int,
  metric_long,
  metric_long_long,
  metric_double,
  metric_string
} metric_type;

static int metric_ticker = 0;
static metric_type types[METRIC_MAX];
static int ints[METRIC_MAX];
static long longs[METRIC_MAX];
static long long long_longs[METRIC_MAX];
static double doubles[METRIC_MAX];
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

void metric_load_long(char *key, long value){
  if (metric_ticker >= METRIC_MAX) {
    printf("ReproBLAS error: too many metrics\n");
    exit(125);
  }
  keys[metric_ticker] = key;
  types[metric_ticker] = metric_long;
  longs[metric_ticker] = value;
  metric_ticker++;
}

void metric_load_long_long(char *key, long long value){
  if (metric_ticker >= METRIC_MAX) {
    printf("ReproBLAS error: too many metrics\n");
    exit(125);
  }
  keys[metric_ticker] = key;
  types[metric_ticker] = metric_long_long;
  long_longs[metric_ticker] = value;
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
  printf("{\n");
  for(i = 0; i < metric_ticker; i++){
    switch(types[i]){
      case metric_int:
        printf("\t\"%s\": %d", keys[i], ints[i]);
        break;
      case metric_long:
        printf("\t\"%s\": %ld", keys[i], longs[i]);
        break;
      case metric_long_long:
        printf("\t\"%s\": %lld", keys[i], long_longs[i]);
        break;
      case metric_double:
        printf("\t\"%s\": %g", keys[i], doubles[i]);
        break;
      case metric_string:
        printf("\t\"%s\": \"%s\"", keys[i], strings[i]);
        break;
    }
    if(i < metric_ticker - 1){
      printf(",\n");
    }else{
      printf("\n");
    }
  }
  printf("}\n");
  metric_ticker = 0;
}
