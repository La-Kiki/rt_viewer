#pragma once

#include <glm/glm.hpp>

#include <cstdlib>
#include <iostream>

namespace rt {
    struct HitRecord;

class Material {
  public:
     virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const = 0;
};

class Lambertian : public Material {
public:
    Lambertian(const glm::vec3 &a) : albedo(a) {}

    virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override {
        auto scatter_direction = rec.normal + random_in_hemisphere(rec.normal); //glm::normalize(rt::random_in_unit_sphere());
        if (rt::near_zero(scatter_direction)) {
            scatter_direction = rec.normal;
        }
            
        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo;

        return true;
    }

public:
    glm::vec3 albedo;
};

class Metal : public Material {
public:
    Metal(const glm::vec3 &a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

    virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override {
        glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()), rec.normal);
        scattered = Ray(rec.p, reflected + (float)fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (glm::dot(scattered.direction(), rec.normal) > 0);
    }

public:
    glm::vec3 albedo;
    double fuzz;
};

class Dielectric : public Material {
public:
    Dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(
        const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered
    ) const override {
        bool normal_in_ray_dir = glm::dot(r_in.direction(), rec.normal) > 0.0;

        double refraction_ratio = normal_in_ray_dir ? ir : (1.0 / ir);
        glm::vec3 normal_dir = normal_in_ray_dir ? -rec.normal : rec.normal;

        glm::vec3 unit_direction = glm::normalize(r_in.direction());
        
        double cos_theta = fmin(glm::dot(-unit_direction, normal_dir), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        glm::vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double()) {
            // Must Reflect
            direction = glm::reflect(unit_direction, normal_dir);
        }
        else {
            // Can Refract
            direction = refract(unit_direction, normal_dir, refraction_ratio);
        }

        attenuation = glm::vec3(1.0, 1.0, 1.0);
        scattered = Ray(rec.p, direction);
        return true;
    }

public:
    double ir; // Index of Refraction

private:
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

}  // namespace rt