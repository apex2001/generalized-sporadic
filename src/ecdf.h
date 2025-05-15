#ifndef ECDF_H
#define ECDF_H

#include <algorithm>
#include <vector>

#include "generate_task_set.h"

// will be removed
// struct Task {
//     int minSeparation_Ti;
//     string criticalityLevel_Li;
//     double wcetVals_Ci[2] = {0, 0}; //[0] = CiL, [1] = CiH
//     double deadline_Di;
//     double deadlineLC_DiL; //DiL = Di for LC tasks and DiL <= Di for HC tasks
// };


//------------ MOD FUNCTION ------------------

double MOD(double t, double Ti) {
    double res = t - floor(t / Ti) * Ti;
    
    return res;
}

//---------- PROPOSITION 1 FUNCTIONS - EDF schedulability in LC-------------

double dbfLi(int time, int Ti, double DiL, double CiL) { //Equation (1)
    double difference = double(time) - double(DiL);
    double divRes = difference / Ti;
    double RHS = (floor(divRes) +  1) * CiL;
    double res = max(0.0, RHS);

    return res;
}

// bool Proposition1(int time, vector<Task> taskSet) {
//     double sum = 0;
//     //int* Ti, int* DiL, int* CiL

//     for (int i = 0; i < taskSet.size(); i++) {
//         sum += dbfLi(time, taskSet[i].minSeparation_Ti, taskSet[i].deadlineLC_DiL, taskSet[i].wcetVals_Ci[0]);
//     }

//     // cout << "Prop1 time: " << time << endl;
//     // cout << "Prop1 sum: " << sum << endl;
//     // cout << endl;
    

//     return sum <= time;
// }

bool Proposition1(TaskSet& taskSet) {
    int t_max = taskSet.get_t_max();
    if (t_max == 0) return false;

    for (int t = 0; t <= t_max; t++) {
        // cout << "i: " << t << endl;
        double sum = 0;
        for (const auto& pair : taskSet.get_task_set()) {
            const Task& task = pair.second;
            sum += dbfLi(t, task.T, task.tight_D, task.C_LO);
        }

        if (sum > t) {
            return false;
        }
    }
    
    return true;
}

double dbfHi_l4(int t1, int t2, double CiH, double Di, int Ti) {
    double RHS = (floor((t2 - t1 - Di) / Ti) + 1) * CiH;

    return max(double(0), RHS);
}

double dbfUNi_l2(int t1, int t2, double CiL, double DiL, double Ti) {
    double res = 0;
    double modRes = MOD(t1, Ti);
    double divRes = t1 / Ti;
    double clause2LHS = (floor(divRes)) * Ti + DiL;

    if (DiL > modRes && clause2LHS <= t2)
        res = min(double(CiL), modRes);

    return res;
}

double dbfUN_l6(int t1, int t2, TaskSet& taskSet) {
    double res = 0;
    double dbfUNi_Sum = 0;
    double max_DiL = 0;
    map<int, Task> taskSetMap = taskSet.get_task_set();

    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::LO || (task.L == Level::HI && ((t2 - t1) <= (task.D - task.tight_D)))) { //or HC tasks and case 1
            dbfUNi_Sum += dbfUNi_l2(t1, t2, task.C_LO,  task.tight_D, task.T);

            if (max_DiL < task.tight_D) {
                max_DiL = task.tight_D;
            }
        }
    }

    res = min(max_DiL, dbfUNi_Sum);

    return res;
}

double dbfL1_l7(int t1, int t2, TaskSet& taskSet) {
    double res = 0;
    double dbfLi_Sum = 0;
    double dbfUN_Res = dbfUN_l6(t1, t2, taskSet);

    map<int, Task> taskSetMap = taskSet.get_task_set();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        // cout << "task T: " << task.T << endl;
        // cout << "task C_LO: " << task.C_LO << endl;
        // cout << task.L << endl;
        if (task.L == Level::LO || (task.L == Level::HI && ((t2 - t1) <= (task.D - task.tight_D)))) {
            dbfLi_Sum += dbfLi(t1, task.T, task.tight_D, task.C_LO);
            // cout << "task T: " << task.T << endl;
        }
    }

    // cout << "dbfLi_Sum " << dbfLi_Sum << endl;
    // cout << "dbfUN_Res " << dbfUN_Res << endl;

    return dbfLi_Sum + dbfUN_Res;
}

double dbfLi_l4(int t1, int t2, double CiL, double Di, int Ti) {
    double RHS = (floor((t2 - Di) / Ti) - floor ((t2 - t1 - Di) / Ti) - 1);

    return max(double(0), RHS) * CiL;
}

double CO_prop2(int time, double CiL, double Di, double DiL, int Ti) {
    double RHS = MOD(time, Ti) - (Di - DiL);

    return min(double(CiL), RHS);
}

double dbfL2_l8(int t1, int t2, TaskSet& taskSet) {
    double res = 0;
    map<int, Task> taskSetMap = taskSet.get_task_set();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        // cout << "task T: " << task.T << endl;
        // cout << "task C_LO: " << task.C_LO << endl;
        // cout << task.L << endl;
        if (task.L == Level::HI && ((t2 - t1) > (task.D - task.tight_D))) {
            double cal = dbfLi_l4(t1, t2, task.C_LO, task.D, task.T)\
                        + task.C_LO\
                        + CO_prop2(t2 - t1, task.C_LO, task.D, task.tight_D, task.T);
            res += cal;
        }
    }

    return res;
}

double dbfL3_l9(int t1, int t2, TaskSet& taskSet) {
    double res = 0;
    map<int, Task> taskSetMap = taskSet.get_task_set();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        // cout << "task T: " << task.T << endl;
        // cout << "task C_LO: " << task.C_LO << endl;
        // cout << task.L << endl;
        if (task.L == Level::HI && ((t2 - t1) > (task.D - task.tight_D))) {
            double cal = dbfLi_l4(t1, t2, task.C_LO, task.D, task.T)\
                        + task.C_LO;
            res += cal;
        }
    }

    return res;
}

double Theorem2_LHS(int t1, int t2, TaskSet& taskSet) {
    double dbfL1_Res = dbfL1_l7(t1, t2, taskSet);
    double dbfL2_Res = dbfL2_l8(t1, t2, taskSet);
    double dbfL3_Res = dbfL3_l9(t1, t2, taskSet);

    double dbfHi_Res = 0; //dbfHi_l4(t1, t2, CiH[0], Di[0], Ti[0]);
    double CO_Res = 0; //CO(t2 - t1, CiL[0], Di[0], DiL[0], Ti[0]);

    map<int, Task> taskSetMap = taskSet.get_task_set();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        // cout << "task T: " << task.T << endl;
        // cout << "task C_LO: " << task.C_LO << endl;
        // cout << task.L << endl;
        if (task.L == Level::HI) {
            if ((t2 - t1) > (task.D - task.tight_D)) { // case 2 only
                double CO_Cal = CO_prop2(t2 - t1, task.C_LO, task.D, task.tight_D, task.T)\
                            + task.C_HI - task.C_LO; //should check if task is HC and case 2
                CO_Res += CO_Cal;
            }

            if ((t2 - t1) > (task.D - task.tight_D)) { // case 2 or 3, confer with prof
                dbfHi_Res += dbfHi_l4(t1, t2, task.C_HI, task.D, task.T);
            }
        }
    }

    double sumDBFs = dbfL1_Res + dbfL2_Res + dbfL3_Res;
    double min_t1 = min(double(t1), sumDBFs);
    double res = min_t1 + dbfHi_Res + CO_Res;

    return res;
}

bool Theorem2(int t1, int t2, TaskSet& taskSet) {
    double res = Theorem2_LHS(t1, t2, taskSet);

    return res <= t2;
} 

int findCandiate(int t1, int t2, TaskSet& taskSet) {
    double Theorem2_LHS_Res = Theorem2_LHS(t1, t2, taskSet);
    double DEM = Theorem2_LHS_Res - t2;
    int res = -1;
    double DIFF = 0, DEC = INFINITY;

    
}

string ECDF(TaskSet& taskSet) {
    return "Success";
}



#endif