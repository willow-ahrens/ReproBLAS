#ifndef __TEST_METRIC_H
#define __TEST_METRIC_H


void metric_load_int(char *key, int value);
void metric_load_double(char *key, double value);
void metric_load_string(char *key, char *value);

void metric_dump();

#endif
