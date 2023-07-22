#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>

#include "car_model.h"
#include "comfun.h"
#include "nmpc1.h"

#include<typeinfo> //查看数据类型

int main()
{
    Eigen::Vector3f states; //小车实时状态
    Eigen::Vector2f real_controls; //真实控制量
    states<<0,0,0;
    real_controls<<1,0.2;

    float ts = 0.2; //仿真离散时间

    Car car(ts,states); //创建小车


    int k=0;  //仿真步数，随仿真时间累加
    float sim_time; //仿真时间

    std::vector<float> y_simdata; //小车运行轨迹
    std::vector<float> x_simdata;

    std::vector<float> t_simdata;
    std::vector<float> v_simdata;
    std::vector<float> omega_simdata; //存储求解出的最优控制序列

    int predict_step=100;  //NMPC预测步长
    NMPC1 nmpc(predict_step,ts); //创建NMPC控制器
 
    std::vector<Eigen::Matrix<float,3,1>> predict_trajectory; //存储NMPC预测轨迹

    Eigen::Vector3f desired_states;  //小车期望状态（控制目标）
    Eigen::Vector2f desired_controls; //小车期望控制量（控制目标）
    Eigen::Vector2f controls;  //最优序列的第一个控制量

    desired_states<<30,30,M_PI/2;
    desired_controls<<0,0;
    controls<<1,0.1;

    while(true)
    {
        //仿真时间递增
        k++;
        sim_time = (k-1)*ts;

        t_simdata.push_back(sim_time); //存储仿真时间

        //NMPC
        nmpc.opti_solution(states,desired_states,desired_controls); //优化求解
        controls = nmpc.get_controls(); //获得最优序列的第一个控制量

        v_simdata.push_back(controls[0]);
        omega_simdata.push_back(controls[1]);

        predict_trajectory = nmpc.get_predict_trajectory();  //获得预测的轨迹

        std::vector<float> predict_x_data; //预测轨迹
        std::vector<float> predict_y_data;

        for(int j=0;j<predict_step+1;j++)
        {
            predict_x_data.push_back(predict_trajectory.at(j)[0]);
            predict_y_data.push_back(predict_trajectory.at(j)[1]);
        }

        // std::cout<<predict_trajectory<<endl;

        // predict_x_data.push_back((casadi::DM) predict_trajectory(0,1));

        //模型状态更新
        car.state_update(controls);
        states = car.get_states();
        real_controls = car.get_controls();

        x_simdata.push_back(states[0]);
        y_simdata.push_back(states[1]);

        //仿真可视化
        if(k%10==0)
        {
            plt::figure(1);
            plt::clf();
            plt::plot(y_simdata,x_simdata,{{"color", "red"}, {"linestyle", "-"}}); 
            plt::plot(predict_y_data,predict_x_data,{{"color", "blue"}, {"linestyle", "--"}});
            //实时绘制小车模型
            Eigen::Vector2f pos;
            float psi; 
            pos<<states(0),states(1);
            psi = states(2);
            std::string color1 = "red";
            std::string color2 = "blue";
            COMFUN::myWAMVplot(pos,psi,color1);
            //绘制小车期望位姿
            Eigen::Vector2f desired_pos;
            float desired_psi;
            desired_pos<< desired_states[0],desired_states[1];
            desired_psi = desired_states[2];
            COMFUN::myWAMVplot(desired_pos,desired_psi,color2);
            float linewid = 1;
            // plt::xlim(-100, 100);
            // plt::ylim(-100, 100);

            plt::figure(2);
            plt::plot(t_simdata,v_simdata,{{"color", "red"}, {"linestyle", "-"}});

            plt::figure(3);
            plt::plot(t_simdata,omega_simdata,{{"color", "red"}, {"linestyle", "-"}});

            plt::pause(0.1);
        }


    }

    //测试
    // casadi::SX mm = casadi::SX::sym("mm",2,2);
    // mm(0,0) = (float)(1);
    // mm(0,1) = (float)(2);
    // mm(1,0) = (float)(3);
    // mm(1,1) = (float)(4);

    // std::cout<< "mm: "<<mm<<endl;

    // casadi::SX xx = casadi::SX::reshape(mm,-1,1);

    // std::cout<<"xx: "<<xx<<endl;

    // std::cout<<"mm: "<<mm<<endl;
    // std::vector<casadi::SX> xx;
    // xx.push_back(mm);
    // std::cout<<"xx: "<<xx<<endl;
    // float a = 10;
    // std::cout<<"数据类型: "<<typeid(xx).name()<<endl;

    // std::cout<< "hello world!"<< std::endl;
}