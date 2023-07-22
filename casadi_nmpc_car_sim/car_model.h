#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <math.h>

class Car
{
private:
    /* data */
    float m_sample_time; //采样时间
    Eigen::Vector3f m_initial_states; //初始状态
    Eigen::Vector3f m_states; //实时状态
    Eigen::Vector2f m_controls; //实时控制

    float v_max; // 最大前向速度【物理约束】
    float omega_max; // 最大转动角速度 【物理约束】
    
    
public:
    Car(float sample_time, Eigen::Vector3f initial_states);
    ~Car();
    //状态更新
    void state_update(Eigen::Vector2f controls);
    //获取状态
    Eigen::Vector3f get_states();
    //获取控制
    Eigen::Vector2f get_controls();
};
