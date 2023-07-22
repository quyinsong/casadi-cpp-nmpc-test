#include "car_model.h"
#include "comfun.h"

Car::Car(float sample_time, Eigen::Vector3f initial_states)
{
    m_sample_time = sample_time;
    m_states = initial_states;
    m_controls<<0,0;

    v_max = 0.6; // 最大前向速度【物理约束】
    omega_max = (float)M_PI/(float)4.0; // 最大转动角速度 【物理约束】
}

Car::~Car()
{
}

void Car::state_update(Eigen::Vector2f controls)
{
    //物理限制
    if(std::abs(controls[0])>v_max)
    {
        if(std::abs(controls[0]>v_max)){
            m_controls(0) = v_max;
        }else{
            m_controls(0) = -v_max;
        }
    }else{
        m_controls(0) = controls[0];
    }

    if(std::abs(controls[1])>omega_max)
    {
        if(std::abs(controls[1]>omega_max)){
            m_controls(1) = omega_max;
        }else{
            m_controls(1) = -omega_max;
        }
    }else{
        m_controls(1) = controls[1];
    }
    Eigen::Vector3f rhs;
    rhs(0) = m_controls[0]*std::cos(m_states[2]);
    rhs(1) = m_controls[0]*std::sin(m_states[2]);
    rhs(2) = m_controls[1];
    m_states = m_states+rhs*m_sample_time;
}

Eigen::Vector3f Car::get_states()
{
    return m_states;
}

Eigen::Vector2f Car::get_controls()
{
    return m_controls;
}
