#ifndef ECDF_H
#define ECDF_H

#include <algorithm>
#include <vector>
#include <cmath>
#include <limits.h>

#include "generate_task_set.h"

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

bool Proposition1(int t, TaskSet& taskSet) {
    // int t_max = taskSet.get_t_max();
    // if (t_max == 0) return false;

    // for (int t = 0; t <= t_max; t++) {
    //     // cout << "i: " << t << endl;
    //     double sum = 0;
    //     for (const auto& pair : taskSet.get_task_set()) {
    //         const Task& task = pair.second;
    //         sum += dbfLi(t, task.T, task.tight_D, task.C_LO);
    //     }

    //     if (sum > t) {
    //         return false;
    //     }
    // }

    
    double sum = 0;
    for (const auto& pair : taskSet.get_task_set()) {
        const Task& task = pair.second;
        sum += dbfLi(t, task.T, task.tight_D, task.C_LO);
    }

    if (sum > t) {
        return false;
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
            if ((t2 - t1) > (task.D - task.tight_D)) {
                if  ((task.D - task.tight_D) < MOD(t2 - t1, task.T) && // case 2 lemma 4
                    (MOD(t2 - t1, task.T) < task.D) &&      
                    ((floor((t2 - t1) / task.T) * task.T + task.D) <= t2)) { // case 2 only

                    double CO_Cal = CO_prop2(t2 - t1, task.C_LO, task.D, task.tight_D, task.T)\
                                + (task.C_HI - task.C_LO); //should check if task is HC and case 2
                    CO_Res += CO_Cal;
                }

                else if (!((task.D - task.tight_D) < MOD(t2 - t1, task.T) && // case 2 lemma 4
                        (MOD(t2 - t1, task.T) < task.D) &&      
                        ((floor((t2 - t1) / task.T) * task.T + task.D) <= t2))  ||
                    
                        !((task.tight_D > MOD(t1, task.T)) && 
                        ((floor(t1 / task.T) / task.T + task.tight_D) <= t2))) { // case 3 lemma 5
                            dbfHi_Res += dbfHi_l4(t1, t2, task.C_HI, task.D, task.T);
                }
            }

            

            // if ((t2 - t1) > (task.D - task.tight_D)) { // case 2 or 3, confer with prof
                
            // }
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

int findCandiate(vector<Task> candidates, int t1, int t2, TaskSet& taskSet) {
    double Theorem2_LHS_Res = Theorem2_LHS(t1, t2, taskSet);
    double DEM = Theorem2_LHS_Res - t2;
    int res = -1;
    double DIFF = 0, DEC = INFINITY;

    for (int i = 0; i < candidates.size(); i++) {
        if (((t2 - t1) > (candidates[i].D - candidates[i].tight_D)) &&   //1st claude of case 2

            ((candidates[i].D - candidates[i].tight_D) < MOD(t2 - t1, candidates[i].T)) && //2nd clause of case 2 lemma 4
            (MOD(t2 - t1, candidates[i].T) < candidates[i].D) &&      


            ((floor((t2 - t1) / candidates[i].T) * candidates[i].T + candidates[i].D) <= t2) &&//2nd clause of case 2 lemma 4

            ((candidates[i].C_HI - candidates[i].C_LO) >= DEM)) {  
        
            if ((MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D)) < DEC) {
                DEC = MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D);
                res = i;
                DIFF = candidates[i].C_HI - candidates[i].C_LO;
            }
            else if ((MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D)) == DEC) {
                if (candidates[i].C_HI - candidates[i].C_LO > DIFF) {
                    res = i;
                    DIFF = candidates[i].C_HI - candidates[i].C_LO;
                }
            }
        }    
    }

    return res;
}

string ECDF(TaskSet& taskSet) {
    int i = -1, tMax = taskSet.get_t_max(), min_diff = INT_MAX;
    vector<Task> candidates;

    map<int, Task> taskSetMap = taskSet.get_task_set();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::HI) {
            candidates.push_back(task);
            min_diff = min(min_diff, (task.D - task.tight_D));
        }
    }
    
    while (true) {
        bool feasible = true;

        for (int t = 0; t <= tMax; t++) {
            if (!Proposition1(t, taskSet)) {
                if (i == -1) {
                    return "Failure";
                }
                
                for (auto& candidate : candidates) {
                    if (candidate.ID == i) {
                        candidate.C_LO += 1;
                        candidates.erase(candidates.begin() + i);
                        break;
                    }
                }
                
                i = -1;
                break;
            }
        }

        for (int t2 = 0; t2 <= tMax; t2++) {
            for (int t1 = 0; t1 < t2 - min_diff; t1++) {
                if (!Theorem2(t1, t2, taskSet)) {
                    if (t1 == 0 && candidates.empty()) {
                        return "Failure";
                    }

                    i = findCandiate(candidates, t1, t2, taskSet);

                    for (auto& candidate : candidates) {
                        if (candidate.ID == i) {
                            candidate.tight_D -= 1;
                            
                            if (candidate.tight_D - 1 < candidate.C_LO) {
                                candidates.erase(candidates.begin() + i);
                            }
                            break;
                        }
                    }

                    feasible = false;
                    break;
                }
            }
        }
        if (feasible) {
            return "Success";
        }
    }

    return "Unknown"; // This line should never be reached
}



#endif