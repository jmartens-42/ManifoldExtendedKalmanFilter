#include <manifoldextendedkalmanfilter/tools.hpp>

Eigen::Matrix<float, 3, 3> crossProdMat(const Eigen::Vector<float, 3>& v){

    Eigen::Matrix<float, 3, 3> m;

    m <<    0, -1*v(2), v(1),
            v(2), 0, -1*v(0),
            -1*v(1), v(0), 0;

    
    return m;

}


Eigen::Quaternionf angularVelocityToRotation(Eigen::Vector<float, 3> w, float dt){

    w = w * (3.14159f/180.0f);
    float norm_over_2 = w.norm()*dt*0.5f;
    float real_part = cosf(norm_over_2);
    auto vector_part = Eigen::Vector<float, 3>(w.normalized()*sinf(norm_over_2));

    // std::cout << "real part: " << real_part << " vector part: " << vector_part.transpose() << " norm_over_2 " << norm_over_2 << "\n";
    return Eigen::Quaternionf(real_part, vector_part.x(), vector_part.y(), vector_part.z()).normalized();

}
