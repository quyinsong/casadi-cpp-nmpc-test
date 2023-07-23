#include "nmpc.h"
#include "asv.h"

NMPC::NMPC(int predict_step, float ts, int integration_step, std::vector<Eigen::Matrix<float,8,1>> trajectory_points)
{
    m_predict_step = predict_step;
    m_sample_time = ts;
    m_integration_step = integration_step;
    m_dt = m_sample_time/m_integration_step;

    m_trajectory_points = trajectory_points;
    m_predict_stage = -1;
}

NMPC::~NMPC()
{
}

//预测模型
casadi::SX NMPC::pridictModel(casadi::SX last_states_init, casadi::SX current_Tl, casadi::SX current_Tr)
{
    //下一时刻状态与上一时刻状态
    casadi::SX next_states = casadi::SX::sym("next_states",6); 
    casadi::SX last_states = casadi::SX::sym("last_states",6); 
    for(int j=0; j<6;j++)
    {
        last_states(j) = last_states_init(j);
    }

    for(int j=0;j<m_integration_step;j++)
    {
        //实时获取无人艇状态
        casadi::SX u = last_states(0);
        casadi::SX v = last_states(1);
        casadi::SX r = last_states(2);
        casadi::SX psi = last_states(5);
        casadi::SX U = sqrt(pow(u,2)+pow(v,2));
        //无人船参数
        casadi::SX Xudot = -11.25, Nvdot = 0.2, Iz = 400, Yvdot = -195.6398;
        casadi::SX Nrdot = -600.2102, xg = 0, Yrdot = 0.2;
        casadi::SX Xu = -50.5870, Xuu = -5.8722;
        casadi::SX Yv = -20*rho*abs(v)*2.1092;
        casadi::SX Nr = -(float)(0.02)*rho*m_pi*U*pow(T,2)*pow(L,2);
        casadi::SX Nv = -(float)(0.06)*rho*m_pi*U*pow(T,2)*L;
        casadi::SX Yr = -6*rho*m_pi*U*pow(T,2)*L;
        casadi::SX Yvv = -599.3130, Yrv = -524.3989, Yvr = -524.3989, Yrr = -1378;
        casadi::SX Nvv = -524.3989, Nrv = -1378, Nvr = -1378, Nrr = -2996;
        casadi::SX m11 = m-Xudot;
        casadi::SX m22 = m-Yvdot;
        casadi::SX m23 = m*xg-Yrdot;
        casadi::SX m32 = m*xg-Nvdot;
        casadi::SX m33 = Iz-Nrdot;
        casadi::SX c13 = -m*v+(Yvdot*v+(Yrdot+Nvdot)*r/2)/200;
        casadi::SX c23 = m*u-Xudot*u;
        casadi::SX c31 = -c13, c32 = -c23;
        casadi::SX d11 = -Xu-Xuu*abs(u);
        casadi::SX d22 = -Yv-Yrv*abs(v)-Yvr*abs(r);
        casadi::SX d23 = -Yr-Yrv*abs(v)-Yrr*abs(r);
        casadi::SX d32 = -Nv-Nvv*abs(v)-Nvr*abs(r);
        casadi::SX d33 = -Nr-Nrv*abs(v)-Nrr*abs(r);
        //系统矩阵描述
        casadi::SX Rbn = casadi::SX::sym("Rbn",3,3);
        Rbn(0) = casadi::SX::cos(psi); Rbn(3) = -casadi::SX::sin(psi); Rbn(6) = 0;
        Rbn(1) = casadi::SX::sin(psi); Rbn(4) = casadi::SX::cos(psi); Rbn(7) = 0;
        Rbn(2) = 0; Rbn(5) = 0; Rbn(8) = 1;
        
        casadi::SX M = casadi::SX::sym("M",3,3);
        M(0) = m11; M(3) = 0; M(6) = 0;
        M(1) = 0; M(4) = m22; M(7) = m23;
        M(2) = 0; M(5) = m32; M(8) = m33;

        casadi::SX Crb = casadi::SX::sym("Crb",3,3);
        Crb(0) = 0; Crb(3) = 0; Crb(6) = c13;
        Crb(1) = 0; Crb(4) = 0; Crb(7) = c23;
        Crb(2) = c31; Crb(5) = c32; Crb(8) = 0;

        casadi::SX Dv = casadi::SX::sym("Dv",3,3);
        Dv(0) = d11; Dv(3) = 0; Dv(6) = 0;
        Dv(1) = 0; Dv(4) = d22; Dv(7) = d23;
        Dv(2) = 0; Dv(5) = d32; Dv(8) = d33;

        casadi::SX nu = casadi::SX::sym("nu",3);
        nu(0) = u; nu(1) = v; nu(2) = r;

        //无人艇状态更新
        casadi::SX BT = casadi::SX::sym("BT",3,2);
        BT(0) = 1;   BT(3) = 1;
        BT(1) = 0;   BT(4) = 0;
        BT(2) = B/2; BT(5) = -B/2;
        casadi::SX tau = casadi::SX::sym("tau",3);
        tau(0) = current_Tl+current_Tr;
        tau(1) = 0;
        tau(2) = (current_Tl-current_Tr)*BT(2);
        casadi::SX nudot = casadi::SX::sym("nudot",3);
        // nudot = casadi::SX::mtimes(casadi::SX::inv(M),
        //         -casadi::SX::mtimes(Crb        // std::cout<<"nudot_size"<<nudot.size()<<endl;


        // std::cout<<"nudot: "<<nudot<<endl;,nu)-casadi::SX::mtimes(Dv,nu)+tau);
        
        nudot(0) = (-c13*r-d11*u+tau(0))/m11;
        nudot(1) = (-c23*r-d22*v-d23*r+tau(1))/m22;
        nudot(2) = (-c31*u-c32*v-d32*v-d33*r+tau(2))/m33;

        casadi::SX etadot = casadi::SX::sym("etadot",3);
        // etadot = casadi::SX::mtimes(Rbn,nu);
        
        etadot(0) = u*casadi::SX::cos(psi)-v*casadi::SX::sin(psi);
        etadot(1) = u*casadi::SX::sin(psi)+v*casadi::SX::cos(psi);
        etadot(2) = r;

        casadi::SX xdot = casadi::SX::sym("xdot",6);
        xdot(0) = nudot(0); xdot(1) = nudot(1); xdot(2) = nudot(2);
        xdot(3) = etadot(0); xdot(4) = etadot(1); xdot(5) = etadot(2);

        next_states = xdot*m_dt+last_states;

        last_states = next_states;
    }
    
    return next_states;
}
//滚动优化求解控制序列
Eigen::Vector2f NMPC::nmpc(Eigen::Matrix<float,6,1> current_states)
{
    //调用次数加1
    m_predict_stage++;
    //上一状态和下一状态声明
    casadi::SX last_states = casadi::SX::sym("last_states",6);
    casadi::SX next_states = casadi::SX::sym("next_states",6);
    // 优化求解的变量
    casadi::SX Tl = casadi::SX::sym("Tl",m_predict_step);
    casadi::SX Tr = casadi::SX::sym("Tr",m_predict_step);
    //代价函数
    casadi::SX cost_fun = casadi::SX::sym("cost_fun",1);
    cost_fun = 0; 
    //每一次预测的误差
    casadi::SX error_states = casadi::SX::sym("error_states",6);
    casadi::SX error_Tl = casadi::SX::sym("error_Tl",1);
    casadi::SX error_Tr = casadi::SX::sym("error_Tr",1);
    //获取预测模型的初始状态
    casadi::SX states =casadi::SX::sym("states",6);
    for(int j=0;j<=5;j++)
    {
        states(j) = current_states[j];
    }
    // std::cout<<last_states<<endl;
    //预测
    for(int k=0;k<m_predict_step;k++)
    {
        //一步预测-----------------------------------------------
        // next_states = pridictModel(last_states, Tl(k), Tr(k));

        for(int j=0;j<m_integration_step;j++)
        {
            //实时获取无人艇状态
            casadi::SX u = states(0);
            casadi::SX v = states(1);
            casadi::SX r = states(2);
            casadi::SX psi = states(5);
            casadi::SX U = sqrt(pow(u,2)+pow(v,2));
            //无人船参数
            casadi::SX Xudot = -11.25, Nvdot = 0.2, Iz = 400, Yvdot = -195.6398;
            casadi::SX Nrdot = -600.2102, xg = 0, Yrdot = 0.2;
            casadi::SX Xu = -50.5870, Xuu = -5.8722;
            casadi::SX Yv = -20*rho*abs(v)*2.1092;
            casadi::SX Nr = -(float)(0.02)*rho*m_pi*U*pow(T,2)*pow(L,2);
            casadi::SX Nv = -(float)(0.06)*rho*m_pi*U*pow(T,2)*L;
            casadi::SX Yr = -6*rho*m_pi*U*pow(T,2)*L;
            casadi::SX Yvv = -599.3130, Yrv = -524.3989, Yvr = -524.3989, Yrr = -1378;
            casadi::SX Nvv = -524.3989, Nrv = -1378, Nvr = -1378, Nrr = -2996;
            casadi::SX m11 = m-Xudot;
            casadi::SX m22 = m-Yvdot;
            casadi::SX m23 = m*xg-Yrdot;
            casadi::SX m32 = m*xg-Nvdot;
            casadi::SX m33 = Iz-Nrdot;
            casadi::SX c13 = -m*v+(Yvdot*v+(Yrdot+Nvdot)*r/2)/200;
            casadi::SX c23 = m*u-Xudot*u;
            casadi::SX c31 = -c13, c32 = -c23;
            casadi::SX d11 = -Xu-Xuu*abs(u);
            casadi::SX d22 = -Yv-Yrv*abs(v)-Yvr*abs(r);
            casadi::SX d23 = -Yr-Yrv*abs(v)-Yrr*abs(r);
            casadi::SX d32 = -Nv-Nvv*abs(v)-Nvr*abs(r);
            casadi::SX d33 = -Nr-Nrv*abs(v)-Nrr*abs(r);
            //系统矩阵描述
            casadi::SX Rbn = casadi::SX::sym("Rbn",3,3);
            Rbn(0) = casadi::SX::cos(psi); Rbn(3) = -casadi::SX::sin(psi); Rbn(6) = 0;
            Rbn(1) = casadi::SX::sin(psi); Rbn(4) = casadi::SX::cos(psi); Rbn(7) = 0;
            Rbn(2) = 0; Rbn(5) = 0; Rbn(8) = 1;
            
            casadi::SX M = casadi::SX::sym("M",3,3);
            M(0) = m11; M(3) = 0; M(6) = 0;
            M(1) = 0; M(4) = m22; M(7) = m23;
            M(2) = 0; M(5) = m32; M(8) = m33;

            casadi::SX Crb = casadi::SX::sym("Crb",3,3);
            Crb(0) = 0; Crb(3) = 0; Crb(6) = c13;
            Crb(1) = 0; Crb(4) = 0; Crb(7) = c23;
            Crb(2) = c31; Crb(5) = c32; Crb(8) = 0;

            casadi::SX Dv = casadi::SX::sym("Dv",3,3);
            Dv(0) = d11; Dv(3) = 0; Dv(6) = 0;
            Dv(1) = 0; Dv(4) = d22; Dv(7) = d23;
            Dv(2) = 0; Dv(5) = d32; Dv(8) = d33;

            casadi::SX nu = casadi::SX::sym("nu",3);
            nu(0) = u; nu(1) = v; nu(2) = r;

            //无人艇状态更新
            casadi::SX BT = casadi::SX::sym("BT",3,2);
            BT(0) = 1;   BT(3) = 1;
            BT(1) = 0;   BT(4) = 0;
            BT(2) = B/2; BT(5) = -B/2;
            casadi::SX tau = casadi::SX::sym("tau",3);
            tau(0) = Tl(k)+Tr(k);
            tau(1) = 0;
            tau(2) = (Tl(k)-Tr(k))*BT(2);
            casadi::SX nudot = casadi::SX::sym("nudot",3);
            
            nudot(0) = (-c13*r-d11*u+tau(0))/m11;
            nudot(1) = (-c23*r-d22*v-d23*r+tau(1))/m22;
            nudot(2) = (-c31*u-c32*v-d32*v-d33*r+tau(2))/m33;

            casadi::SX etadot = casadi::SX::sym("etadot",3);
            // etadot = casadi::SX::mtimes(Rbn,nu);
            
            etadot(0) = u*casadi::SX::cos(psi)-v*casadi::SX::sin(psi);
            etadot(1) = u*casadi::SX::sin(psi)+v*casadi::SX::cos(psi);
            etadot(2) = r;

            casadi::SX xdot = casadi::SX::sym("xdot",6);
            xdot(0) = nudot(0); xdot(1) = nudot(1); xdot(2) = nudot(2);
            xdot(3) = etadot(0); xdot(4) = etadot(1); xdot(5) = etadot(2);

            states = xdot*m_dt+states;
        }
    //--------------------------------------------------------------------------------
        // std::cout<<next_states<<endl;
        // //先测试定点控制
        // casadi::SX desired_point = casadi::SX::sym("desired_point",6);
        // desired_point(0) = 1; desired_point(1) = 0; desired_point(2) = 0.1;
        // desired_point(3) = 20; desired_point(4) = 20; desired_point(5) = 0;
        
        //测试轨迹点跟踪
        casadi::SX current_desired_trajectory = casadi::SX::sym("current_desired_trajectory",6);
        casadi::SX current_desired_inputs = casadi::SX::sym("current_desired_inputs",2);

        int trajectory_points_index;
        trajectory_points_index = m_predict_stage+k;

        std::cout<<"trajectory_points_index: "<<trajectory_points_index<<endl;

        if(trajectory_points_index<m_trajectory_points.size()-1)
        {
            for(int j=0;j<8;j++)
            {
                if(j<6){
                    current_desired_trajectory(j) = m_trajectory_points.at(trajectory_points_index)(j);
                }else{
                    current_desired_inputs(j-6) = m_trajectory_points.at(trajectory_points_index)(j);
                }
            }
        }
        else
        {
            for(int j=0;j<8;j++)
            {
                if(j<6){
                    current_desired_trajectory(j) = m_trajectory_points.at(m_trajectory_points.size()-1)(j);
                }else{
                    current_desired_inputs(j-6) = m_trajectory_points.at(m_trajectory_points.size()-1)(j);
                }
            }
        }
        
        // std::cout<< "m_trajectory_points_size: "<< m_trajectory_points.size()<<endl; 
        //计算误差
        error_states = states;
        error_Tl = Tl(k);
        error_Tr = Tr(k);
        //计算代价函数
        cost_fun = cost_fun+
                   1000*casadi::SX::dot(error_states(0),error_states(0))+
                   1000*casadi::SX::dot(error_states(1),error_states(1))+
                   1000*casadi::SX::dot(error_states(2),error_states(2))+
                   1000*casadi::SX::dot(error_states(3),error_states(3))+
                   1000*casadi::SX::dot(error_states(4),error_states(4))+
                   1000*casadi::SX::dot(error_states(5),error_states(5))+
                   0*casadi::SX::dot(error_Tl,error_Tl)+
                   0*casadi::SX::dot(error_Tr,error_Tr);

        // cost_fun = cost_fun+
        //            100*error_Tl*error_Tl+
        //            100*error_Tr*error_Tr;
        //更新状态 
        // last_states = next_states;

    }

    //约束
    casadi::SX inputs = casadi::SX::vertcat({Tl,Tr}); //待求解变量
    casadi::SX constraints = 0;

    // std::cout<<"inputs: "<< inputs<<endl;
    // std::cout<<"constraints: "<< error_states<<endl;

    // Initial guess and bounds for the optimization variables
    std::vector<float> inputs_initial_value;
    std::vector<float> inputs_min, inputs_max;
    
    for(int j=0; j<2*m_predict_step;j++)
    {
        inputs_initial_value.push_back(100);
        inputs_max.push_back(200.0);
        inputs_min.push_back(-200.0);
    }

    // std::cout<<"inputs_max: "<<inputs_max<<endl;
    // std::cout<<"inputs_min: "<<inputs_min<<endl;


    // Nonlinear bounds
    std::vector<float> constraints_max = { (float)casadi::inf, (float)casadi::inf, (float)casadi::inf,
                        (float) casadi::inf, (float)casadi::inf, (float)casadi::inf };
    std::vector<float> constraints_min = { -(float)casadi::inf, -(float)casadi::inf, 
                        -(float)casadi::inf, -(float)casadi::inf, -(float)casadi::inf, -(float)casadi::inf };

    // std::cout<<"constraints_max: "<<constraints_max<<endl;
    // std::cout<<"constraints_min: "<<constraints_min<<endl;

    // NLP
    casadi::SXDict nlp;                 // NLP declaration
    nlp["x"]= inputs;                   // decision vars
    nlp["f"] = cost_fun;                // objective
    // nlp["g"] = constraints;             // constraints

    // Create NLP solver and buffers
    casadi::Function solver = nlpsol("solver", "ipopt", nlp);
    std::map<std::string, casadi::DM> arg, res;

    // Solve the NLP
    arg["lbx"] = inputs_min;
    arg["ubx"] = inputs_max;
    // arg["lbg"] = -10;
    // arg["ubg"] = 10;
    arg["x0"] = inputs_initial_value;
    res = solver(arg);

    // 获取代价函数和控制序列
    std::vector<float> res_control_all(res.at("x"));

    std::vector<float> res_control_Tl, res_control_Tr;

    res_control_Tl.assign(res_control_all.begin(), res_control_all.begin() + m_predict_step);
    res_control_Tr.assign(res_control_all.begin() + m_predict_step, res_control_all.begin() + 2 * m_predict_step);

    // 采用求解得到的控制序列的第一组作为当前控制量
    Eigen::Vector2f control_input_now;
    control_input_now << res_control_Tl.front(), res_control_Tr.front();
    // control_input_now<< (float)120,(float)100;

    // std::cout<<"error_states: "<<error_states<<endl;
    // std::cout<<"next_states: "<< next_states<<endl;
    std::cout<<"cost_fun: "<< cost_fun<<endl;

    std::cout << "objective: " << res.at("f") << endl;
    std::cout << "solution: " << res.at("x") << endl;

    std::cout<<"res_control_Tl: "<< res_control_Tl<<endl;
    std::cout<<"res_control_Tr: "<< res_control_Tr<<endl;

    std::cout << "control_input_now: " << control_input_now[0]<<" "<< control_input_now[1] << endl;

    return control_input_now;


}