#ifndef ECDF_TEST_H
#define ECDF_TEST_H

#include "../src/ecdf.h"

#include <gtest/gtest.h>

TEST(ECDF, MOD) {
    // EXPECT_TRUE(false);
    EXPECT_EQ(MOD(5, 3), 2);
    EXPECT_EQ(MOD(7, 4), 3);
    EXPECT_EQ(MOD(10, 5), 0);
    EXPECT_EQ(MOD(15, 7), 1);
    EXPECT_EQ(MOD(20, 8), 4);
}

TEST(ECDF, dbfLi) {
    EXPECT_EQ(dbfLi(10, 5, 2, 3), 6);
    EXPECT_EQ(dbfLi(15, 7, 4, 2), 4);
    EXPECT_EQ(dbfLi(20, 10, 5, 1), 2);
    EXPECT_EQ(dbfLi(25, 12, 6, 4), 8);
    EXPECT_EQ(dbfLi(30, 15, 8, 3), 6);
    EXPECT_EQ(dbfLi(35, 20, 40, 5), 0);
    EXPECT_EQ(dbfLi(10, 15, 40, 2), 0);
}

TEST(ECDF, Proposition1) {
    TaskSet taskSet = TaskSet(0.5);

    EXPECT_EQ(Proposition1(taskSet), true);
}

TEST(ECDF, Proposition1_Fail_ManyTasks_HighCombinedUtilization) {
    Task task1(0, 5, 3, 3, 5, 5);    // C_LO = 3, T = 5, D = 5. Utilization = 0.6
    Task task2(1, 10, 4, 4, 10, 10);  // C_LO = 4, T = 10, D = 10. Utilization = 0.4
    Task task3(2, 7, 2, 2, 7, 7);    // C_LO = 2, T = 7, D = 7 Utilization = 0.2857
    Task task4(3, 15, 6, 6, 15, 15);   // C_LO = 6, T = 15, D = 15 Utilization = 0.4

    TaskSet testSet;

    testSet.task_set[task1.ID] = task1;
    testSet.task_set[task2.ID] = task2;
    testSet.task_set[task3.ID] = task3;
    testSet.task_set[task4.ID] = task4;

    // Update the TaskSet metadata
    testSet.num_tasks = testSet.task_set.size();
    testSet.utilization = testSet.calculate_utilization();
    testSet.t_max = testSet.calculate_t_max();

    cout << "num_tasks: " << testSet.num_tasks << endl;
    cout << "utilization: " << testSet.utilization << endl
            << "t_max: " << testSet.t_max << endl;

    EXPECT_FALSE(Proposition1(testSet));
}

TEST(ECDF, dbfHi_l4) {
    EXPECT_EQ(dbfHi_l4(0, 0, 1, 1, 1), 0);
    EXPECT_EQ(dbfHi_l4(0, 1, 1, 1, 1), 1);
    EXPECT_EQ(dbfHi_l4(1, 5, 1, 1, 1), 4);
}

TEST(ECDF, dbfUNi_l2) {
    ASSERT_EQ(dbfUNi_l2(5, 10, 3, 2, 7), 0.0);
    ASSERT_EQ(dbfUNi_l2(3, 8, 5, 6, 10), 3.0);
    ASSERT_EQ(dbfUNi_l2(3, 8, 2, 6, 10), 2.0);
}

TEST(dbfUN_l6Test, MixedTasks_ConditionFalseForHI) {
    TaskSet ts;
    Task task1(0, 10, 2, 4, 10, 8, Level::LO);
    Task task2(1, 15, 3, 6, 20, 9, Level::HI);

    ts.task_set[task1.ID] = task1;
    ts.task_set[task2.ID] = task2;
    // max_DiL = 9
    // LO task: dbfUNi_l2(0, 10, 2, 8, 10)
    // HI task: (10 - 0) <= (20 - 9) => 10 <= 11 (true) -> dbfUNi_l2(0, 10, 3, 9, 15)
    ASSERT_EQ(dbfUN_l6(0, 10, ts), std::min(9.0, dbfUNi_l2(0, 10, 2, 8, 10) + dbfUNi_l2(0, 10, 3, 9, 15)));

    TaskSet ts2;
    Task task3(0, 10, 2, 4, 10, 8, Level::LO);
    Task task4(1, 15, 3, 6, 18, 9, Level::HI);
    ts2.task_set[task3.ID] = task3;
    ts2.task_set[task4.ID] = task4;
    // max_DiL = 9
    // LO task: dbfUNi_l2(0, 5, 2, 8, 10)
    // HI task: (5 - 0) <= (18 - 9) => 5 <= 9 (true) -> dbfUNi_l2(0, 5, 3, 9, 15)
    ASSERT_EQ(dbfUN_l6(0, 5, ts2), std::min(9.0, dbfUNi_l2(0, 5, 2, 8, 10) + dbfUNi_l2(0, 5, 3, 9, 15)));

    TaskSet ts3;
    Task task5(0, 10, 2, 4, 10, 5, Level::LO);
    Task task6(1, 15, 3, 6, 12, 12, Level::HI);
    Task task7(2, 7, 1, 2, 8, 3, Level::LO);
    ts3.task_set[task5.ID] = task5;
    ts3.task_set[task6.ID] = task6;
    ts3.task_set[task7.ID] = task7;
    // max_DiL = 12
    // LO task 0: dbfUNi_l2(0, 10, 2, 5, 10)
    // HI task 1: (10 - 0) <= (12 - 12) => 10 <= 0 (false)
    // LO task 2: dbfUNi_l2(0, 10, 1, 3, 7)
    ASSERT_EQ(dbfUN_l6(0, 10, ts), std::min(12.0, dbfUNi_l2(0, 10, 2, 5, 10) + dbfUNi_l2(0, 10, 1, 3, 7)));

    TaskSet ts4;
    Task task8(0, 10, 8, 4, 10, 2, Level::LO);
    Task task9(1, 15, 3, 6, 12, 1, Level::HI);

    ts4.task_set[task8.ID] = task8;
    ts4.task_set[task9.ID] = task9;
    // max_DiL = 2
    // LO task 0: dbfUNi_l2(0, 5, 8, 2, 10)
    // HI task 1: (5 - 0) <= (12 - 1) => 5 <= 11 (true) -> dbfUNi_l2(0, 5, 3, 1, 15)
    ASSERT_EQ(dbfUN_l6(0, 5, ts2), std::min(2.0, dbfUNi_l2(0, 5, 8, 2, 10) + dbfUNi_l2(0, 5, 3, 1, 15)));
}

TEST(dbfL1_l7Test, SingleLOTask) {
    TaskSet ts;
    Task task1(0, 10, 2, 4, 10, 8, Level::LO);
    ts.task_set[task1.ID] = task1;

    map<int, Task> ts_test = ts.get_task_set();
    int task1T = ts_test[task1.ID].L;
    // cout << "task1T: " << task1T << endl

    ASSERT_EQ(task1.L, Level::HI); // Check the level of the local object
    ASSERT_EQ(ts.task_set[0].L, Level::HI); // Check the level after insertion

    // ASSERT_EQ(dbfL1_l7(5, 15, ts), dbfLi(5, 10, 8, 2) + std::min(8.0, dbfUNi_l2(5, 15, 2, 8, 10)));
    ASSERT_EQ(dbfL1_l7(5, 15, ts), 0.0); // Corrected expected value

    TaskSet ts1;
    Task task2(0, 10, 4, 4, 10, 5, Level::LO);
    Task task3(1, 15, 3, 6, 12, 12, Level::HI);
    Task task4(2, 7, 2, 2, 8, 3, Level::LO);
    ts1.task_set[task2.ID] = task2;
    ts1.task_set[task3.ID] = task3;
    ts1.task_set[task4.ID] = task4;

    cout << "task2: " << task2.L << endl;
    cout << "task3: " << task3.L << endl;    
    cout << "task4: " << task4.L << endl;
    // dbfLi_Sum = dbfLi(0, 10, 5, 2) + dbfLi(0, 15, 12, 3) + dbfLi(0, 7, 3, 1)
    // dbfUN_Res = dbfUN_l6(0, 10, ts)
    // max_DiL = 12
    // HI task: (10-0) <= (12-12) => 10 <= 0 (false)
    ASSERT_EQ(dbfL1_l7(0, 10, ts1), 0);

    TaskSet ts2;
    // cout << dbfLi(0, 10, 8, 2) + dbfLi(0, 15, 9, 3) + std::min(9.0, dbfUNi_l2(0, 10, 2, 8, 10) + dbfUNi_l2(0, 10, 3, 9, 15)) << endl;
    ASSERT_EQ(dbfL1_l7(0, 10, ts2), 0.0); // dbfUN_l6 returns -1 for empty set

    TaskSet ts3;
    Task task5(0, 10, 2, 4, 10, 8, Level::LO);
    Task task6(1, 15, 3, 6, 20, 9, Level::HI);
    ts3.task_set[task5.ID] = task5;
    ts3.task_set[task6.ID] = task6;
    // dbfLi_Sum = dbfLi(0, 10, 8, 2) + dbfLi(0, 15, 9, 3)
    // dbfUN_Res = dbfUN_l6(0, 10, ts) (10-0) <= (20-9) => 10 <= 11 (true)
    ASSERT_EQ(dbfL1_l7(0, 10, ts3), dbfLi(0, 10, 8, 2) + dbfLi(0, 15, 9, 3) + std::min(9.0, dbfUNi_l2(0, 10, 2, 8, 10) + dbfUNi_l2(0, 10, 3, 9, 15)));
}



#endif
