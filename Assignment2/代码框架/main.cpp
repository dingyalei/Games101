#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_rotation(Eigen::Vector3f axis, float angel){
    Eigen::Matrix4f rotation = Eigen::Matrix4f::Identity();

    float x = axis.x();
    float y = axis.y();
    float z = axis.z();
    float sin = std::sin(angel);
    float cos = std::cos(angel);
    rotation << x*x*(1-cos)+cos, x*y*(1-cos)+z*sin, x*z*(1-cos)-y*sin, 0,
                x*y*(1-cos)-z*sin, y*y*(1-cos)+cos, y*z*(1-cos)+x*sin, 0,
                x*z*(1-cos)+y*sin, y*z*(1-cos)-x*sin, z*z*(1-cos)+cos,0,
                0,0,0,1;
    return rotation;


}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    float angle = std::sin(rotation_angle/180.0 * acos(-1));
    // Eigen::Matrix4f rotation_matrix;
    // model << std::cos(angle), std::sin(-angle), 0, 0,
    //          std::sin(angle), std::cos(angle),  0, 0,
    //          0,               0,                1, 0,
    //          0,               0,                0, 1;


    return get_rotation(Eigen::Vector3f{0,0,1},angle);
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    
    // projection <<std::tan(MY_PI/(eye_fov/2))/aspect_ratio, 0, 0, 0,
    //           0, std::tan(MY_PI/(-eye_fov/2)), 0, 0,
    //           0, 0, (zNear+zFar)/(zNear-zFar), -2*zNear*zFar/(zNear-zFar),
    //           0, 0, 1, 0;

    
    // float l = aspect_ratio * zNear * std::tan(eye_fov/2);
    // float r = -aspect_ratio * zNear * std::tan(eye_fov/2);
    // float t = zNear * tan(eye_fov/2);
    // float b = -zNear * tan(eye_fov/2);
    // float n = -zNear;
    // float f = -zFar;

    // Eigen::Matrix4f translate;
    // Eigen::Matrix4f scale;
    // translate<< 1, 0, 0, -(l+r)/2.0f,
    //             0, 1, 0, -(b+t)/2.0f,
    //             0, 0, 1, -(f+n)/2.0f,
    //             0, 0, 0, 1;
    // scale<< 2.0f/(l-r), 0, 0, 0,
    //         0, 2.0f/(t-b), 0, 0,
    //         0, 0, 2.0f/(n-f), 0,
    //         0, 0, 0, 1.0;
    // proj<< n, 0, 0, 0,
    //         0, n, 0, 0,
    //         0, 0, n+f, -n*f,
    //         0, 0, 1, 0;
    // // projection = scale * translate * proj * projection;
    // projection = scale * projection;
    // float eye_fov_rad = eye_fov / 180.0f * acos(-1);
    // float t = zNear * tan(eye_fov_rad/2.0f);
    // float r = t * aspect_ratio;
    // float b = -t;
    // float l = -r;
    // float n = -zNear;
    // float f = -zFar;
    // // frustum -> cubic
    // Eigen::Matrix4f M_p2o;
    // M_p2o << n, 0, 0, 0,
    //          0, n, 0, 0,
    //          0, 0, n+f, -n*f,
    //          0, 0, 1, 0;
    // // orthographic projection
    // Eigen::Matrix4f M_o_shift = Eigen::Matrix4f::Identity();
    // M_o_shift(0, 3) = -(r+l)/2.0f;
    // M_o_shift(1, 3) = -(t+b)/2.0f;
    // M_o_shift(2, 3) = -(n+f)/2.0f;
    // Eigen::Matrix4f M_o_scale = Eigen::Matrix4f::Identity();
    // M_o_scale(0, 0) = 2.0f / (r-l);
    // M_o_scale(1, 1) = 2.0f / (t-b);
    // M_o_scale(2, 2) = 2.0f / (n-f);
    // // squash all transformations
    // projection = M_o_scale * M_o_shift * M_p2o * projection;
    // std::clog << "projection" << std::endl << projection << std::endl;

    // float top = -tan((eye_fov/2.0f) * abs(zNear));
    // float right = top * aspect_ratio;

    // projection << zNear/right,0,0,0,
    //               0,zNear/top,0,0,
    //               0,0,(zNear+zFar)/(zNear-zFar),(2*zNear*zFar)/(zFar-zNear),
    //               0,0,1,0;
     Eigen::Matrix4f P2O = Eigen::Matrix4f::Identity();
    P2O << zNear, 0, 0, 0,
        0, zNear, 0, 0,
        0, 0, zNear + zFar, (-1) * zFar * zNear,
        0, 0, 1, 0;
    float halfEyeAngelRadian = eye_fov / 2.0 / 180.0 * MY_PI;
    float t = zNear * std::tan(halfEyeAngelRadian); //top y轴的最高点
    float r = t * aspect_ratio;                     //right x轴的最大值
    float l = (-1) * r;                             //left x轴最小值
    float b = (-1) * t;                             //bottom y轴的最大值
    //进行一定的缩放使之成为一个标准的长度为2的正方体
    Eigen::Matrix4f ortho1 = Eigen::Matrix4f::Identity();
    ortho1 << 2 / (r - l), 0, 0, 0,
        0, 2 / (t - b), 0, 0,
        0, 0, 2 / (zNear - zFar), 0,
        0, 0, 0, 1;
    // 把一个长方体的中心移动到原点
    Eigen::Matrix4f ortho2 = Eigen::Matrix4f::Identity();
    ortho2 << 1, 0, 0, (-1) * (r + l) / 2,
        0, 1, 0, (-1) * (t + b) / 2,
        0, 0, 1, (-1) * (zNear + zFar) / 2,
        0, 0, 0, 1;
    Eigen::Matrix4f Matrix_ortho = ortho1 * ortho2;
    projection = Matrix_ortho * P2O;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
