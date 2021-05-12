#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <cstdlib>
#include <iostream>

#include <vector>

class Lambertian;
namespace rt {



struct RTContext {
    int width = 500;
    int height = 500;
    std::vector<glm::vec4> image;
    bool freeze = false;
    int current_frame = 0;
    int current_line = 0;
    int max_frames = 1000;
    int max_bounces = 1;
    float epsilon = 2e-4f;
    glm::mat4 view = glm::mat4(1.0f);
    glm::vec3 ground_color = glm::vec3(1.0f, 1.0f, 1.0f);
    glm::vec3 sky_color = glm::vec3(0.5f, 0.7f, 1.0f);
    bool show_normals = true;
    // Add more settings and parameters here
    std::shared_ptr<Lambertian> material_ground = std::make_shared<Lambertian>(glm::vec3(0.8, 0.8, 0.0));
    /*glm::vec3 material_center = std::make_shared<Lambertian>(glm::vec3(0.7, 0.3, 0.3));
    glm::vec3 material_left = std::make_shared<Metal>(glm::vec3(0.8, 0.8, 0.8));
    glm::vec3 material_right = std::make_shared<Metal>(glm::vec3(0.8, 0.6, 0.2));*/
    // ...
};

void setupScene(RTContext &rtx, const char *mesh_filename);
void updateImage(RTContext &rtx);
void resetImage(RTContext &rtx);
void resetAccumulation(RTContext &rtx);

glm::vec3 random_in_unit_sphere();
bool near_zero(glm::vec3 e);

}  // namespace rt
