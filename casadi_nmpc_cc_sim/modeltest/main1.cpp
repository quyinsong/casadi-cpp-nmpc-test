#include "asv.h"
#include "comfun.h"
#include <casadi/casadi.hpp>
#include "nmpc.h"
#include "single_trajectory.h"
#include "nmpc1.h"

using namespace casadi;

//无人艇仿真模型测试
int main()
{
    float ts = 0.02; //采样步长
    float total_time = 150; //仿真总时长
    int Ns = total_time/ts;
    int integration_step = 1;
    Eigen::Matrix<float,6,1> x_init;
    x_init << 0,0,0,0,0,m_pi/2;
    WAMV14 asv1(ts,integration_step, x_init); //创建无人艇对象
    Eigen::Matrix<float,3,1> tau_w;
    Eigen::Matrix<float,2,1> Thrust_c;
    tau_w<<0,0,0;
    Thrust_c<<120,100;
    std::vector<float> x_simdata, y_simdata, psi_simdata; //存储无人艇运动状态
    Eigen::Matrix<float,2,1> pos; //无人艇位置向量
    float psi; //无人艇艏向角
    //存储时间
    std::vector<float> t_sim_data;
    std::vector<float> Tu_sim_data;
    std::vector<float> Tr_sim_data;

    //生成轨迹点
    Eigen::Matrix<float,6,1> x_init_ttg;
    x_init_ttg << 0,0,0,10,10,m_pi/2;
    SingleTTG singlettg(ts, total_time, integration_step,  x_init_ttg);
    std::vector<Eigen::Matrix<float,8,1>> trajectory_points;
    trajectory_points = singlettg.TTG();
    // std::cout<< trajectory_points.at(800)<<endl;
    std::vector<float> trajectory_points_x;
    std::vector<float> trajectory_points_y;
    std::vector<float> trajectory_points_psi;
    std::cout<<"trajectory_points.size(): " << trajectory_points.size()<<endl;
    for(int k=0;k<trajectory_points.size();k++)
    {
      trajectory_points_x.push_back(trajectory_points.at(k)(3));
      trajectory_points_y.push_back(trajectory_points.at(k)(4));
      trajectory_points_psi.push_back(trajectory_points.at(k)(5));
    }
    Eigen::Matrix<float,2,1> trajectory_pos; //轨迹点位置向量
    float trajectory_psi; //轨迹点位置向量

    //设定期望位置
    Eigen::Matrix<float,6,1> desired_states;
    desired_states<< 0,0,0,50,50,0;

    //创建NMPC对象
    int predict_step = 20;
    // NMPC m_nmpc(predict_step, ts, integration_step, trajectory_points);
    // NMPC1 m_nmpc1(predict_step, ts, desired_states);
    NMPC1 m_nmpc1(predict_step, ts, trajectory_points);

    //优化控制输出
    Eigen::Vector2f tau;
    tau<<0,0;

    //预测轨迹
    std::vector<Eigen::Matrix<float,6,1>> predict_trajectoy;


    int k=0;
    while(true)
    {   

        k++;

        float t = ts*(k-1);
        //存储仿真时间
        t_sim_data.push_back(t);
        // cout<< t<<endl;
        //获取无人艇状态
        Eigen::Matrix<float,6,1> x;
        x = asv1.getstate();
        //获取无人艇推力
        Eigen::Matrix<float,2,1> thrust;
        thrust = asv1.getthrust();
        //存储无人艇运动状态
        x_simdata.push_back(x(3));
        y_simdata.push_back(x(4));
        psi_simdata.push_back(x(5));

        //调用nmpc求解最优控制量
        // Thrust_c = m_nmpc.nmpc(x);
        m_nmpc1.opti_solution(x);

        //获得最优控制量
        tau = m_nmpc1.get_controls();

        //存储推力
        Tu_sim_data.push_back(tau[0]);
        Tr_sim_data.push_back(tau[1]);

        //获得预测轨迹
        std::vector<float> predict_x;
        std::vector<float> predict_y;
        predict_trajectoy = m_nmpc1.get_predict_trajectory();
        for(int j=0;j<predict_step+1;j++)
        {
            predict_x.push_back(predict_trajectoy.at(j)[3]);
            predict_y.push_back(predict_trajectoy.at(j)[4]);
        }

        //控制力坐标变换
        Eigen::Matrix<float,2,2> BT;
        BT << 1,1,1.83/2,-1.83/2;
        Thrust_c = BT.inverse()*tau;

        //无人艇状态更新
        asv1.state_update(Thrust_c,tau_w);

        //仿真可视化
        if(k%100==0)
        {
            plt::figure(1);
            plt::clf();
            plt::plot(y_simdata,x_simdata,{{"color", "red"}, {"linestyle", "-"}});  
            plt::plot(trajectory_points_y,trajectory_points_x,{{"color", "blue"}, {"linestyle", "-"}});  
            pos<<x(3),x(4);
            psi = x(5);
            if(k<trajectory_points.size()-1)
            {
                trajectory_pos<<trajectory_points_x.at(k-1), trajectory_points_y.at(k-1);
                trajectory_psi = trajectory_points_psi.at(k-1);
            }else{
                trajectory_pos<<trajectory_points_x.at(trajectory_points_x.size()-1),
                                trajectory_points_y.at(trajectory_points_y.size()-1);
                trajectory_psi = trajectory_points_psi.at(trajectory_points_psi.size()-1);

            }
            std::string color1 = "red";
            std::string color2 = "blue";
            COMFUN::myWAMVplot(pos,psi,color1);

            //绘制预测轨迹
            plt::plot(predict_y,predict_x,{{"color", "blue"}, {"linestyle", "--"}});  

            //绘制期望位置
            Eigen::Vector2f desired_pos;
            float desired_psi;
            desired_pos<<desired_states[3],desired_states[4];
            desired_psi = desired_states[5];
            // COMFUN::myWAMVplot(desired_pos,desired_psi,color2);
            COMFUN::myWAMVplot(trajectory_pos,trajectory_psi,color2);

            float linewid = 1;
            // plt::xlim(-100, 100);
            // plt::ylim(-100, 100);
            plt::pause(0.1);

            plt::figure(2);
            plt::plot(t_sim_data,Tu_sim_data,{{"color", "red"}, {"linestyle", "-"}});
            plt::plot(t_sim_data,Tr_sim_data,{{"color", "blue"}, {"linestyle", "-"}});
            plt::pause(0.1);
        }


    }

    // casadi::SX XX = casadi::SX::sym("XX",6);

    // casadi::SX YY = XX(casadi::Slice(1,6,1));

    // std::cout<<"YY: "<< YY<<endl;

    // //测试casadi矩阵乘法
    // casadi::SX xxx = casadi::SX::sym("xxx",3,3);
    // casadi::SX yyy = casadi::SX::sym("yyy",3);
    // casadi::SX zzz = casadi::SX::mtimes(xxx,yyy);
    // std::cout<< "zzz: "<<zzz<< endl;

//     //以下代码测试casadi的NLP求解
//     cout << "casadi_test" << endl;

//     casadi::SX Rbn = casadi::SX::sym("Rbn",3,3);
//     Rbn(0) = 1; Rbn(1) = 2;
    

// // This is another way to define a nonlinear solver. Opti is new
// /*
//   *    min  x1*x4*(x1 + x2 + x3) + x3
//   *    s.t. x1*x2*x3*x4 >=25
//             x1^2 + x2^2 + x3^2 + x4^2 = 40
//             1 <= x1, x2, x3, x4 <= 5
// */

// //   Optimization variables
//   SX x = SX::sym("x", 4);
//   std::cout << "x:" << x << std::endl;
  
//   // Objective
//   SX f = x(0)*x(3)*(x(0) + x(1) + x(2)) + x(2)+2;
//   std::cout << "f:" << f << std::endl;

//   // Constraints
//   SX g = vertcat(x(0)*x(1)*x(2)*x(3), pow(x(0),2) + pow(x(1),2) + pow(x(2),2) + pow(x(3),2));
//   std::cout << "g:" << g << std::endl;

//   // Initial guess and bounds for the optimization variables
//   vector<double> x0 = { 0.0, 0.0, 0.0, 0.0 };
//   vector<double> lbx = { 1, 1, 1, 1 };
//   vector<double> ubx = {5, 5, 5, 5 };

//   // Nonlinear bounds
//   vector<double> lbg = { 25, 40 };
//   vector<double> ubg = { inf, 40 };

//   // NLP
//   SXDict nlp = { { "x", x }, { "f", f }, { "g", g } };

//   // Create NLP solver and buffers
//   Function solver = nlpsol("solver", "ipopt", nlp);
//   std::map<std::string, DM> arg, res;

//   // Solve the NLP
//   arg["lbx"] = lbx;
//   arg["ubx"] = ubx;
//   arg["lbg"] = lbg;
//   arg["ubg"] = ubg;
//   arg["x0"] = x0;
//   res = solver(arg);

//   // Print the solution
//   cout << "--------------------------------" << endl;
//   // std::cout << res << std::endl;
    // cout << "objective: " << res.at("f") << endl;
    // cout << "solution: " << res.at("x") << endl;

//   return 0;
}
