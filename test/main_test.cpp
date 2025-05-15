#include "gtest/gtest.h"

// #include "test_generate_tasks.h"
// #include "test_generate_task_set.h"
// #include "test_dbfs.h"
// #include "test_schedulability.h"
// #include "test_search_algorithm.h"
// #include "test_ekberg.h"
// #include "test_amc.h"
// #include "test_edf_vd.h"
#include "test_ecdf.h"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}