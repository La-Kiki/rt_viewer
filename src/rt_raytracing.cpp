#include "rt_raytracing.h"
#include "rt_ray.h"
#include "rt_hitable.h"
#include "rt_sphere.h"
#include "rt_triangle.h"
#include "rt_box.h"
#include "material.h"

#include "cg_utils2.h"  // Used for OBJ-mesh loading
#include <stdlib.h>     // Needed for drand48()


namespace rt {

// Store scene (world) in a global variable for convenience
struct Scene {
    Sphere ground;
    std::vector<Sphere> spheres;
    std::vector<Box> boxes;
    std::vector<Sphere> bounding_spheres;
    std::vector<Triangle> mesh;
    Box mesh_bbox;
    std::vector<Box> boundingBoxes;
    std::vector<std::shared_ptr<Material>> material_ptr;//vector with material pointers
} g_scene;

bool hit_world(RTContext &rtx, const Ray &r, float t_min, float t_max, HitRecord &rec)
{
    HitRecord temp_rec;
    bool hit_anything = false;
    float closest_so_far = t_max;

    if (g_scene.ground.hit(r, t_min, closest_so_far, temp_rec)) {
        hit_anything = true;
        closest_so_far = temp_rec.t;
        rec = temp_rec;
    }
    if (rtx.showSpheres) {
        for (int i = 0; i < g_scene.spheres.size(); ++i) {
            if (g_scene.spheres[i].hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
    }
   if (rtx.showBoxes) {
        for (int i = 0; i < g_scene.boxes.size(); ++i) {
            if (g_scene.boxes[i].hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
    }

    if (rtx.showMesh) {
        /* FOR BOUNDING SPHERES
        for (int i = 0; i < g_scene.bounding_spheres.size(); ++i) {
            if (g_scene.bounding_spheres[i].hit(r, t_min, closest_so_far, temp_rec)) {
            */
        for (int i = 0; i < g_scene.boundingBoxes.size(); ++i) {
            if (g_scene.boundingBoxes[i].hit(r, t_min, closest_so_far, temp_rec)) {
                for (int j = 0; j < g_scene.mesh.size(); ++j) {
                    if (g_scene.mesh[j].hit(r, t_min, closest_so_far, temp_rec)) {
                        hit_anything = true;

                        closest_so_far = temp_rec.t;
                        rec = temp_rec;
                    }
                }
            }
        }
        
    }
    return hit_anything;
}

double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max - min) * random_double();
}

glm::vec3 random(double min, double max) {
    return glm::vec3(random_double(min, max), random_double(min, max), random_double(min, max));
}

// Taken from chapter 8 in the "Ray Tracing in a Weekend" book
glm::vec3 random_in_unit_sphere() {
    while (true) {
        auto p = random(-1, 1);
        if (glm::dot(p, p) >= 1) continue;
        return p;
    }
}

glm::vec3 random_in_hemisphere(const glm::vec3& normal) {
    glm::vec3 in_unit_sphere = random_in_unit_sphere();
    if (glm::dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

bool near_zero(glm::vec3 e) {
    // Return true if the vector is close to zero in all dimensions.
    const auto s = 1e-8;
    return (fabs(e.x) < s) && (fabs(e.y) < s) && (fabs(e.z) < s);
}


glm::vec3 refract(const glm::vec3& uv, const glm::vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(glm::dot(-uv, n), 1.0f);
    glm::vec3 r_out_perp = (float)etai_over_etat * (uv + (float)cos_theta * n);
    glm::vec3 r_out_parallel = -sqrt((float)fabs(1.0f - glm::dot(r_out_perp, r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}


// This function should be called recursively (inside the function) for
// bouncing rays when you compute the lighting for materials, like this
//
// if (hit_world(...)) {
//     ...
//     return color(rtx, r_bounce, max_bounces - 1);
// }
//
// See Chapter 7 (+8.2) in the "Ray Tracing in a Weekend" book
glm::vec3 color(RTContext &rtx, const Ray &r, int max_bounces)
{
    if (max_bounces < 0) return glm::vec3(0.0f);

    HitRecord rec;
    if (hit_world(rtx, r, 0.001f, 9999.0f, rec)) {
        rec.normal = glm::normalize(rec.normal);  // Always normalise before use!
        if (rtx.show_normals) {
            return rec.normal * 0.5f + 0.5f; 
        }
        

        Ray scattered;
        glm::vec3 attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation * color(rtx, scattered, max_bounces - 1);
        }

        return glm::vec3(0, 0, 0);
    }

    // If no hit, return sky color
    glm::vec3 unit_direction = glm::normalize(r.direction());
    float t = 0.5f * (unit_direction.y + 1.0f);
    return (1.0f - t) * rtx.ground_color + t * rtx.sky_color;
}


// MODIFY THIS FUNCTION!
void setupScene(RTContext &rtx, const char *filename)
{
    auto materialGround = std::make_shared<Lambertian>(glm::vec3(rtx.groundSphereColor));
    auto materialLambert = std::make_shared<Lambertian>(rtx.lambertianColor);
    auto materialGlass = std::make_shared<Dielectric>(rtx.dielectricRefraction);
    auto materialYellowMetal = std::make_shared<Metal>(glm::vec3(0.8, 0.6, 0.2), rtx.metalFuzz);
    auto materialGreyMetal = std::make_shared<Metal>(glm::vec3(0.77, 0.78, 0.82), rtx.metalFuzz);
    
    //Adds the material pointer to a vector with material pointers
    // First glass element will be used for solid spheres, the second for hollow ones
    g_scene.material_ptr.push_back(materialLambert);
    g_scene.material_ptr.push_back(materialYellowMetal);
    g_scene.material_ptr.push_back(materialGreyMetal);
    g_scene.material_ptr.push_back(materialGlass);
    g_scene.material_ptr.push_back(materialGlass);
    
    

    g_scene.ground = Sphere(glm::vec3(0.0f, -1000.5f, 0.0f), 1000.0f, materialGround);
    g_scene.spheres = {
        Sphere(glm::vec3(2.0f, -0.3f, 0.0f), 0.2f, materialLambert),
        Sphere(glm::vec3(-1.0f, -0.1f, -1.0f), -0.4f, materialGlass),
        Sphere(glm::vec3(-1.0f, -0.2f, -0.6f), 0.3f, materialGlass),
        Sphere(glm::vec3(-1.5f, 0.0f, 0.5f), 0.5f, materialYellowMetal),
    };
    //Enum corresponds to material_ptr index for each sphere in g_scene - 
    // the material used when creating sphere
    rtx.sphereMaterials.push_back(rtx.LAMBERTIAN);
    rtx.sphereMaterials.push_back(rtx.HOLLOWGLASS);
    rtx.sphereMaterials.push_back(rtx.GLASS);
    rtx.sphereMaterials.push_back(rtx.YELLOWMETAL);
    

    
    g_scene.boxes = {
        Box(glm::vec3(0.0f, -0.5f, 1.0f), glm::vec3(0.25f), materialYellowMetal),
        Box(glm::vec3(1.0f, -0.25f, -0.5f), glm::vec3(0.25f), materialGlass),
        Box(glm::vec3(-1.0f, -0.5f, 1.0f), glm::vec3(0.25f), materialGreyMetal),
    };
    rtx.boxMaterials.push_back(rtx.YELLOWMETAL);
    rtx.boxMaterials.push_back(rtx.GLASS);
    rtx.boxMaterials.push_back(rtx.GREYMETAL);
    

        cg::OBJMesh mesh;
        cg::objMeshLoad(mesh, filename);
        g_scene.mesh.clear();

        glm::vec3 maxVert;
        glm::vec3 minVert;

        for (int i = 0; i < mesh.indices.size(); i += 3) {
            int i0 = mesh.indices[i + 0];
            int i1 = mesh.indices[i + 1];
            int i2 = mesh.indices[i + 2];

            glm::vec3 v0 = mesh.vertices[i0] + glm::vec3(0.0f, 0.135f, 0.0f);
            glm::vec3 v1 = mesh.vertices[i1] + glm::vec3(0.0f, 0.135f, 0.0f);
            glm::vec3 v2 = mesh.vertices[i2] + glm::vec3(0.0f, 0.135f, 0.0f);
            g_scene.mesh.push_back(Triangle(v0, v1, v2, materialGreyMetal));

            maxVert = glm::max(glm::max(glm::max(maxVert, v0), v1), v2);
            minVert = glm::min(glm::min(glm::min(minVert, v0), v1), v2);
        }
        rtx.meshMaterial = rtx.GREYMETAL;

        glm::vec3 meshPos = minVert + maxVert * (1 / 2.0f);
        //float meshRadius = glm::length(maxVert - minVert) / 2;   //For bounding circle
        glm::vec3 meshRadius = maxVert - minVert * (1/ 2.0f); //For bounding box

        // Bounding sphere roughly same size and placement as mesh below 
       // Multiple bounding spheres can be added for the same mesh if needed
        //g_scene.bounding_spheres = { Sphere(meshPos, meshRadius)};
        g_scene.boundingBoxes = { Box(meshPos, meshRadius) };
        
}

void updateSpheres(RTContext& rtx) {
    for (int i = 0; i < g_scene.spheres.size(); ++i) {

        // Updates the sphere's current material 
        g_scene.spheres[i].mat_ptr = g_scene.material_ptr[rtx.sphereMaterials[i]];

        // Sets the radius of a sphere to positive 
        g_scene.spheres[i].radius = abs(g_scene.spheres[i].radius);

        // If the material is a dielectric shell the radius is set to negative
        if (rtx.sphereMaterials[i] == rtx.HOLLOWGLASS && (g_scene.spheres[i].radius > 0.0))
        {
            g_scene.spheres[i].radius *= -1;
        }
    }
}

void updateBoxes(RTContext& rtx) {
    for (int i = 0; i < g_scene.boxes.size(); ++i) {

        // Updates the box's current material 
        g_scene.boxes[i].mat_ptr = g_scene.material_ptr[rtx.boxMaterials[i]];

        // Sets the radius of a box to positive 
        g_scene.boxes[i].radius = glm::abs(g_scene.boxes[i].radius);

        bool positiveRadius = glm::all(glm::greaterThan(g_scene.boxes[i].radius, glm::vec3(0)));
        // If the material is a dielectric shell the radius is set to negative
        if (rtx.boxMaterials[i] == rtx.HOLLOWGLASS && positiveRadius)
        {
            g_scene.boxes[i].radius *= -1;
        }
    }
}

void updateMesh(RTContext& rtx, const char* filename) {

    cg::OBJMesh mesh;
    cg::objMeshLoad(mesh, filename);
    g_scene.mesh.clear();

    for (int i = 0; i < mesh.indices.size(); i += 3) {
        int i0 = mesh.indices[i + 0];
        int i1 = mesh.indices[i + 1];
        int i2 = mesh.indices[i + 2];

        glm::vec3 v0 = mesh.vertices[i0] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v1 = mesh.vertices[i1] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v2 = mesh.vertices[i2] + glm::vec3(0.0f, 0.135f, 0.0f);
        g_scene.mesh.push_back(Triangle(v0, v1, v2, g_scene.material_ptr[rtx.meshMaterial]));

    }
}


void updateMaterialPtrs(RTContext &rtx) {
    g_scene.material_ptr[rtx.LAMBERTIAN].reset( new Lambertian(rtx.lambertianColor) );

    g_scene.material_ptr[rtx.YELLOWMETAL].reset( new Metal(glm::vec3(0.8, 0.6, 0.2), rtx.metalFuzz) );
    g_scene.material_ptr[rtx.GREYMETAL].reset( new Metal(glm::vec3(0.77, 0.78, 0.82), rtx.metalFuzz) );

    g_scene.material_ptr[rtx.GLASS].reset( new Dielectric(rtx.dielectricRefraction) );
    g_scene.material_ptr[rtx.HOLLOWGLASS].reset( new Dielectric(rtx.dielectricRefraction) );

    // Choose correct constructor for ground material pointer
    if (rtx.groundMaterial == rtx.LAMBERTIAN) {
        g_scene.ground.mat_ptr.reset(new Lambertian(rtx.groundSphereColor));
    }
    else if (rtx.groundMaterial == rtx.YELLOWMETAL) {
        g_scene.ground.mat_ptr.reset( new Metal(rtx.metalColor1, rtx.metalFuzz) );
    }
    else if (rtx.groundMaterial == rtx.GREYMETAL) {
        g_scene.ground.mat_ptr.reset( new Metal(rtx.metalColor2, rtx.metalFuzz) );
    }
    else if (rtx.groundMaterial == rtx.GLASS || rtx.groundMaterial == rtx.HOLLOWGLASS) {
        g_scene.ground.mat_ptr.reset( new Dielectric(rtx.dielectricRefraction) );
    }
    
}


// MODIFY THIS FUNCTION!
void updateLine(RTContext &rtx, int y, const char* filename)
{
    int nx = rtx.width;
    int ny = rtx.height;
    float aspect = float(nx) / float(ny);
    glm::vec3 lower_left_corner(-1.0f * aspect, -1.0f, -1.0f);
    glm::vec3 horizontal(2.0f * aspect, 0.0f, 0.0f);
    glm::vec3 vertical(0.0f, 2.0f, 0.0f);
    glm::vec3 origin(0.0f, 0.0f, 0.0f);
    glm::mat4 world_from_view = glm::inverse(rtx.view);

    updateMaterialPtrs(rtx);

    //Updates the material of the spheres
    updateSpheres(rtx);

    //Updates the material of the boxes
    updateBoxes(rtx);

    if (rtx.showMesh) {
        if (g_scene.mesh.front().mat_ptr != g_scene.material_ptr[rtx.meshMaterial]) {
            updateMesh(rtx, filename);
        }
    }
        

    // You can try parallelising this loop by uncommenting this line:
    #pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < nx; ++x) {

        float u = (float(x) + 0.5) / (float(nx));
        float v = (float(y) + 0.5) / (float(ny));

        if (rtx.antiAliasingOn)
        {
            u = (float(x) + random_double()) / (float(nx));
            v = (float(y) + random_double()) / (float(ny));
        }
        
        Ray r(origin, lower_left_corner + u * horizontal + v * vertical);
        r.A = glm::vec3(world_from_view * glm::vec4(r.A, 1.0f));
        r.B = glm::vec3(world_from_view * glm::vec4(r.B, 0.0f));

        // Note: in the RTOW book, they have an inner loop for the number of
        // samples per pixel. Here, you do not need this loop, because we want
        // some interactivity and accumulate samples over multiple frames
        // instead (until the camera moves or the rendering is reset).

        if (rtx.current_frame <= 0) {
            // Here we make the first frame blend with the old image,
            // to smoothen the transition when resetting the accumulation
            glm::vec4 old = rtx.image[y * nx + x];
            rtx.image[y * nx + x] = glm::clamp(old / glm::max(1.0f, old.a), 0.0f, 1.0f);
        }
        glm::vec3 c = color(rtx, r, rtx.max_bounces);

        rtx.image[y * nx + x] += glm::vec4(c, 1.0f);
    }
}

void updateImage(RTContext &rtx, const char* filename)
{
    if (rtx.freeze) return;                    // Skip update
    rtx.image.resize(rtx.width * rtx.height);  // Just in case...

    updateLine(rtx, rtx.current_line % rtx.height, filename);

    if (rtx.current_frame < rtx.max_frames) {
        rtx.current_line += 1;
        if (rtx.current_line >= rtx.height) {
            rtx.current_frame += 1;
            rtx.current_line = rtx.current_line % rtx.height;
        }
    }
}

void resetImage(RTContext &rtx)
{
    rtx.image.clear();
    rtx.image.resize(rtx.width * rtx.height);
    rtx.current_frame = 0;
    rtx.current_line = 0;
    rtx.freeze = false;
}

void resetAccumulation(RTContext &rtx)
{
    rtx.current_frame = -1;
}

}  // namespace rt
