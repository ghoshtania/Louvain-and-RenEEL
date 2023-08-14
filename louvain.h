#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "body.h"
//void gcopy(struct insert I,struct insert *Ia);
//void freeG(struct insert I);
struct part Louvain(Graph inputG, struct inserting I);
