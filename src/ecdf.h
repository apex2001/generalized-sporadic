#ifndef ECDF_H
#define ECDF_H

#include <algorithm>
#include <vector>
#include <cmath>
#include <limits.h>

#include "generate_task_set.h"

//------------ MOD FUNCTION ------------------

int MOD(int t, int Ti) {
    double divRes = double(t) / double(Ti);
    int res = t - floor(divRes) * Ti;
    
    return res;
}

//---------- PROPOSITION 1 FUNCTIONS - EDF schedulability in LC-------------

int dbfLi(int time, int Ti, int DiL, int CiL) { //Equation (1)
    int difference = time - DiL;
    double divRes = double(difference) / double(Ti);
    int RHS = (floor(divRes) +  1) * CiL;
    int res = max(0, RHS);

    return res;
}

bool Proposition1(int t, TaskSet& taskSet) {
    int sum = 0;

    // for (int i = 0; i < candidates.size(); i++) {
    //     Task task = candidates[i];
    //     sum += dbfLi(t, task.T, task.tight_D, task.C_LO);
    // }
    for (const auto& pair : taskSet.get_task_set_ref()) {
        const Task& task = pair.second;
        sum += dbfLi(t, task.T, task.tight_D, task.C_LO);
    }

    if (sum > t) {
        return false;
    }
    
    
    return true;
}

int dbfHi_l4(int t1, int t2, int CiH, int Di, int Ti) {
    int interval = t2 - t1;
    int difference = interval - Di;
    double divRes = double(difference) / double(Ti);
    int RHS = (floor(divRes) + 1) * CiH;

    return max(0, RHS);
}

int dbfUNi_l2(int t1, int t2, int CiL, int DiL, int Ti) {
    int res = 0;
    int modRes = MOD(t1, Ti);
    double divRes = double(t1) / double(Ti);
    int clause2LHS = (floor(divRes)) * Ti + DiL;

    if (DiL > modRes && clause2LHS <= t2)
        res = min(CiL, modRes);

    return res;
}

bool isCase1(int t1, int t2, Task task) {
    return (task.L == Level::HI && ((t2 - t1) <= (task.D - task.tight_D)));
}

bool isCase2(int t1, int t2, Task task) {
    return (task.L == Level::HI && ((t2 - t1) > (task.D - task.tight_D)) &&
            ((task.D - task.tight_D) < MOD(t2 - t1, task.T)) && // case 2 lemma 4 pass
            (MOD(t2 - t1, task.T) < task.D) &&      
            ((floor((t2 - t1) / task.T) * task.T + task.D) <= t2));
}

//TODO: check if gemini version vs this version is correct
bool isCase3(int t1, int t2, Task task) {
    return (task.L == Level::HI && ((t2 - t1) > (task.D - task.tight_D)) &&
            !isCase2(t1, t2, task)  &&
            
                !((task.tight_D > MOD(t1, task.T)) && // case 3 lemma 2 fail
                ((floor(t1 / task.T) / task.T + task.tight_D) <= t2)));
}


int dbfUN_l6(int t1, int t2, TaskSet& taskSet) {
    int res = 0;
    int dbfUNi_Sum = 0;
    int max_DiL = 0;
    auto& taskSetMap = taskSet.get_task_set_ref();

    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::LO || isCase1(t1, t2, task)) { //or HC tasks and case 1
            dbfUNi_Sum += dbfUNi_l2(t1, t2, task.C_LO,  task.tight_D, task.T);

            if (max_DiL < task.tight_D) {
                max_DiL = task.tight_D;
            }
        }
    }

    res = min(max_DiL, dbfUNi_Sum);

    return res;
}

int dbfL1_l7(int t1, int t2, TaskSet& taskSet) {
    int res = 0;
    int dbfLi_Sum = 0;
    int dbfUN_Res = dbfUN_l6(t1, t2, taskSet);

    auto& taskSetMap = taskSet.get_task_set_ref();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::LO || isCase1(t1, t2, task)) { //or HC tasks and case 1
            dbfLi_Sum += dbfLi(t1, task.T, task.tight_D, task.C_LO);
        }
    }

    res = dbfLi_Sum + dbfUN_Res;

    return res;
}

int dbfLi_l4(int t1, int t2, int CiL, int Di, int Ti) {
    int RHS = (floor((t2 - Di) / Ti) - floor ((t2 - t1 - Di) / Ti) - 1);

    return max(0, RHS) * CiL;
}

int CO_prop2(int time, int CiL, int Di, int DiL, int Ti) {
    int RHS = MOD(time, Ti) - (Di - DiL);

    return min(CiL, RHS);
}

int dbfL2_l8(int t1, int t2, TaskSet& taskSet) {
    int res = 0;
    auto& taskSetMap = taskSet.get_task_set_ref();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (isCase2(t1, t2, task)) {
            int cal = dbfLi_l4(t1, t2, task.C_LO, task.D, task.T)\
                        + task.C_LO\
                        + CO_prop2(t2 - t1, task.C_LO, task.D, task.tight_D, task.T);
            res += cal;
        }
    }

    return res;
}

int dbfL3_l9(int t1, int t2, TaskSet& taskSet) {
    int res = 0;
    auto& taskSetMap = taskSet.get_task_set_ref();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (isCase3(t1, t2, task)) {
            int cal = dbfLi_l4(t1, t2, task.C_LO, task.D, task.T)
                        + task.C_LO;
            res += cal;
        }
    }

    return res;
}

int Theorem2_LHS(int t1, int t2, TaskSet& taskSet) {
    int dbfL1_Res = dbfL1_l7(t1, t2, taskSet);
    int dbfL2_Res = dbfL2_l8(t1, t2, taskSet);
    int dbfL3_Res = dbfL3_l9(t1, t2, taskSet);

    int dbfHi_Res = 0; //dbfHi_l4(t1, t2, CiH[0], Di[0], Ti[0]);
    int CO_Res = 0; //CO(t2 - t1, CiL[0], Di[0], DiL[0], Ti[0]);

    auto& taskSetMap = taskSet.get_task_set_ref();
    
    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::HI) {
            if ((t2 - t1) > (task.D - task.tight_D)) {
                if  (isCase2(t1, t2, task)) { // case 2 lemma 4

                    int CO_Cal = CO_prop2(t2 - t1, task.C_LO, task.D, task.tight_D, task.T)
                                + (task.C_HI - task.C_LO); 
                    CO_Res += CO_Cal;
                    dbfHi_Res += dbfHi_l4(t1, t2, task.C_HI, task.D, task.T);
                }

                else if (isCase3(t1, t2, task)) { // case 3 lemma 2
                    dbfHi_Res += dbfHi_l4(t1, t2, task.C_HI, task.D, task.T);
                }
            }
        }
    }

    int sumDBFs = dbfL1_Res + dbfL2_Res + dbfL3_Res;
    int min_t1 = min(t1, sumDBFs);
    int res = min_t1 + dbfHi_Res + CO_Res;

    return res;
}

bool Theorem2(int t1, int t2, TaskSet& taskSet) {
    int res = Theorem2_LHS(t1, t2, taskSet);

    return res <= t2;
} 

pair<int, int> findCandidate(vector<Task> candidates, int t1, int t2, TaskSet& taskSet) {
    double Theorem2_LHS_Res = Theorem2_LHS(t1, t2, taskSet);
    double DEM = Theorem2_LHS_Res - t2;
    int pos_in_candidates = -1, ID = -1;
    double DIFF = 0, DEC = INFINITY;

    for (int i = 0; i < candidates.size(); i++) {
        if (isCase2(t1, t2, candidates[i]) &&

            ((candidates[i].C_HI - candidates[i].C_LO) >= DEM)) {  
        
            if ((MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D)) < DEC) {
                DEC = MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D);
                pos_in_candidates = i;
                ID = candidates[i].ID;
                DIFF = candidates[i].C_HI - candidates[i].C_LO;
            }
            else if ((MOD(t2 - t1, candidates[i].T) - (candidates[i].D - candidates[i].tight_D)) == DEC) {
                if (candidates[i].C_HI - candidates[i].C_LO > DIFF) {
                    pos_in_candidates = i;
                    ID = candidates[i].ID;
                    DIFF = candidates[i].C_HI - candidates[i].C_LO;
                }
            }
        }    
    }

    return {pos_in_candidates, ID};
} 

void removeCandidateByIndex(std::vector<Task>& candidates, int index) {
    if (index >= 0 && index < candidates.size()) {
        candidates.erase(candidates.begin() + index);
    }
}

void modifyTaskByID(TaskSet& taskSet, int taskId, int newTightD) { 
    std::map<int, Task>& tasksMap = taskSet.get_task_set_ref(); 
    
    auto it = tasksMap.find(taskId);
    if (it != tasksMap.end()) {
        it->second.tight_D = newTightD; 
    } else {
        std::cerr << "Warning: Task with ID " << taskId << " not found in modifyTaskByID." << std::endl;
    }
}


string ECDF(TaskSet& taskSet) {
    int tMax = taskSet.get_t_max(), min_diff = INT_MAX;
    vector<Task> candidates;
    const auto& taskSetMap = taskSet.get_task_set_ref();

    pair<int, int> i_test = {-1, -1};

    for (const auto& pair : taskSetMap) {
        const Task& task = pair.second;
        
        if (task.L == Level::HI) {
            candidates.push_back(task);
            min_diff = min(min_diff, (task.D - task.tight_D));
        }
    }

    // for debugging purposes
    int whileCount = 0;
    // cout << "tMax: " << tMax << "\n";

    while (true) {
        bool feasible = true;
        whileCount++;

        for (int t = 0; t <= tMax; t++) {
            if (!Proposition1(t, taskSet)) {              
                if (i_test.second == -1) { return "Failure"; }
                
                std::map<int, Task>& tasksMap = taskSet.get_task_set_ref();
                auto it_task_in_set = tasksMap.find(i_test.second);

                if (it_task_in_set != tasksMap.end()) {
                    it_task_in_set->second.tight_D += 1;
                } else {
                    std::cerr << "Critical Error: P1 trying to modify non-existent task ID: " << i_test.second << std::endl;
                    return "Error";
                }

                
                for (size_t k = 0; k < candidates.size(); ++k) {
                    if (candidates[k].ID == i_test.second) {
                        removeCandidateByIndex(candidates, k);
                        break; 
                    }
                }
                
                i_test = {-1, -1};
                break;
            }
        }

        bool hc_check_failed_this_iteration = false;
        for (int t2 = 0; t2 <= tMax; t2++) {
            for (int t1 = 0; t1 < t2 - min_diff; t1++) {
                if (!Theorem2(t1, t2, taskSet)) {
                    if (t1 == 0 || candidates.empty()) {
                        return "Failure";
                    }

                    i_test = findCandidate(candidates, t1, t2, taskSet);

                    if (i_test.first >= 0 && i_test.second >= 0) {
                        candidates[i_test.first].tight_D -= 1;
                        modifyTaskByID(taskSet, i_test.second, candidates[i_test.first].tight_D);

                        if (candidates[i_test.first].tight_D < candidates[i_test.first].C_LO) { 
                            removeCandidateByIndex(candidates, i_test.first);
                        }
                    }

                    // candidates[i_test.first].tight_D -= 1;
                    // modifyTaskByID(taskSet, i_test.second, candidates[i_test.first].tight_D);

                    // if (candidates[i_test.first].tight_D < candidates[i_test.first].C_LO) { 
                    //     removeCandidateByIndex(candidates, i_test.first);
                    // }
                    
                    // hc_check_failed_this_iteration = true;
                    // feasible = false;
                    // break;
                }
            }

            if (hc_check_failed_this_iteration) {
                break; 
            }
        }

        if (feasible == true) {
            return "Success";
        }
    }

    return "Unknown"; // This line should never be reached
}

#endif