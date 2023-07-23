#ifndef NMPC_H_
#define NMPC_H_

#include <casadi/casadi.hpp>
#include <iostream>
#include "comfun.h"

class NMPC
{
private:
    /* data */
    int m_predict_step; //一次采样间隔预测的步数
    float m_sample_time; //采样时间
    int m_integration_step; //积分步长
    float m_dt; //单次积分时间
    std::vector<Eigen::Matrix<float,8,1>> m_trajectory_points; 
    int m_predict_stage; //存储第几次调用nmpc

    float Thrustdot_max; //推力速度限幅
    float Thrust_max; //推力限幅
    float Ths; //无人艇推力积分常数
    float L = 4.29;  // length
    float T = 0.127; // draft
    float Bhull = 0.37; // 浮筒beam
    float B = 1.83; // 浮筒中心间距
    float rho = pow(10,3); // 水的密度
    float m = 150; // mass
    float Cd = 1.1; // 方形系数
    float BOA = 2.2; //width
public:
    NMPC(int predict_step, float ts, int integration_step, std::vector<Eigen::Matrix<float,8,1>> trajectory_points);
    ~NMPC();
    //滚动优化求解控制序列
    Eigen::Vector2f nmpc(Eigen::Matrix<float,6,1> current_states);

private:
    //预测模型
    casadi::SX pridictModel(casadi::SX last_states_init, casadi::SX current_Tl, casadi::SX current_Tr);

};



#endif