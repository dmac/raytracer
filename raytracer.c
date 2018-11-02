#include <limits.h>
#include <math.h>
#include <stdbool.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define WIDTH 1024
#define HEIGHT 1024

#define NUM_LIGHTS 1
#define NUM_OBJECTS 2

typedef struct Color { unsigned char r, g, b, a; } Color;

typedef struct Vec3 { float x, y, z; } Vec3;

typedef struct Camera {
    Vec3  pos;
    Vec3  orn;
    float fov;
} Camera;

typedef struct Sphere {
    Vec3  pos;
    Color color;
    float radius;
    float specular;
    float lambert;
    float ambient;
} Sphere;

typedef struct Scene {
    Camera camera;
    Vec3   lights[NUM_LIGHTS];
    Sphere objects[NUM_OBJECTS];
} Scene;

typedef struct Ray {
    Vec3 pos;
    Vec3 orn;
} Ray;

double vec3_dot(Vec3 u, Vec3 v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vec3 vec3_cross(Vec3 u, Vec3 v) {
    return (Vec3){
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x,
    };
}

Vec3 vec3_add(Vec3 u, Vec3 v) {
    return (Vec3){
        u.x + v.x,
        u.y + v.y,
        u.z + v.z,
    };
}

Vec3 vec3_sub(Vec3 u, Vec3 v) {
    return (Vec3){
        u.x - v.x,
        u.y - v.y,
        u.z - v.z,
    };
}

Vec3 vec3_unit(Vec3 v) {
    float mag = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return (Vec3){
        v.x / mag,
        v.y / mag,
        v.z / mag,
    };
}

Vec3 vec3_scale(Vec3 v, double s) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
}

Vec3 vec3_reflect(Vec3 v, Vec3 n) {
    return vec3_sub(v, vec3_scale(n, 2 * vec3_dot(v, n)));
}

double deg2rad(double deg) {
    return M_PI * deg / 180;
}

Vec3 sphere_normal(Sphere *sphere, Vec3 point) {
    return vec3_unit(vec3_sub(point, sphere->pos));
}

bool sphere_intersect(Sphere *sphere, Ray ray, /* out */ double *dist) {
    Vec3   eye_to_center = vec3_sub(sphere->pos, ray.pos);
    double v             = vec3_dot(eye_to_center, ray.orn);
    double eo_dot        = vec3_dot(eye_to_center, eye_to_center);
    double discriminant  = (sphere->radius * sphere->radius) - eo_dot + (v * v);
    if (discriminant < 0) {
        return false;
    }
    *dist = v - sqrt(discriminant);
    return true;
}

Sphere *scene_intersect(Scene *scene, Ray ray, /* out */ double *dist) {
    Sphere *closest = NULL;
    *dist           = MAXFLOAT;
    for (int i = 0; i < NUM_OBJECTS; i++) {
        Sphere *sphere     = &scene->objects[i];
        double  d          = 0;
        bool    intersects = sphere_intersect(sphere, ray, &d);
        if (intersects && d < *dist) {
            closest = sphere;
            *dist   = d;
        }
    }
    return closest;
}

bool scene_light_visible(Scene *scene, Vec3 light, Vec3 point) {
    Ray ray = {
        .pos = point,
        .orn = vec3_unit(vec3_sub(point, light)),
    };
    double  dist;
    Sphere *sphere = scene_intersect(scene, ray, &dist);
    if (sphere == NULL) {
        return false;
    }
    return dist > -0.005;
}

Color scene_trace(Scene *scene, Ray ray, int depth);
Color scene_surface(Scene *scene, Ray ray, Sphere *sphere, Vec3 point_at_time, Vec3 normal, int depth);

Color scene_trace(Scene *scene, Ray ray, int depth) {
    if (depth > 3) {
        return (Color){0, 0, 0, 0};
    };
    double  dist   = 0;
    Sphere *sphere = scene_intersect(scene, ray, &dist);
    if (sphere == NULL) {
        return (Color){255, 255, 255, 255};
    }
    Vec3 point_at_time = vec3_add(ray.pos, vec3_scale(ray.orn, dist));
    return scene_surface(scene, ray, sphere, point_at_time, sphere_normal(sphere, point_at_time), depth);
}

Color scene_surface(Scene *scene, Ray ray, Sphere *sphere, Vec3 point_at_time, Vec3 normal, int depth) {
    Vec3   c       = {0, 0, 0};
    double lambert = 0;
    if (sphere->lambert > 0) {
        for (int i = 0; i < NUM_LIGHTS; i++) {
            Vec3 light = scene->lights[i];
            if (!scene_light_visible(scene, light, point_at_time)) {
                continue;
            }
            double contribution = vec3_dot(vec3_unit(vec3_sub(light, point_at_time)), normal);
            if (contribution > 0) {
                lambert += contribution;
            }
        }
    }
    if (sphere->specular > 0) {
        Ray reflected_ray = {
            .pos = point_at_time,
            .orn = vec3_reflect(ray.orn, normal),
        };
        Color reflected_color = scene_trace(scene, reflected_ray, depth + 1);
        if (reflected_color.a > 0) {
            Vec3 color_vec = {reflected_color.r, reflected_color.g, reflected_color.b};
            c              = vec3_add(c, vec3_scale(color_vec, sphere->specular));
        }
    }
    if (lambert > 1) {
        lambert = 1;
    }
    Vec3 color_sphere = {sphere->color.r, sphere->color.g, sphere->color.b};
    Vec3 color_ret    = vec3_add(c, vec3_scale(color_sphere, lambert * sphere->lambert));
    color_ret         = vec3_add(color_ret, vec3_scale(color_sphere, sphere->ambient));
    return (Color){color_ret.x, color_ret.y, color_ret.z, 255};
}

int main(void) {
    Scene scene  = {0};
    scene.camera = (Camera){
        .pos = {0, 1.8, 10},
        .orn = {0, 3, 0},
        .fov = 45,
    };
    scene.lights[0]  = (Vec3){-30, -10, 20};
    scene.objects[0] = (Sphere){
        .pos      = {0, 3.5, -3},
        .color    = {155, 200, 155, 255},
        .radius   = 3,
        .specular = 0.2,
        .lambert  = 0.7,
        .ambient  = 0.1,
    };
    scene.objects[1] = (Sphere){
        .pos   = {-3, 2, 0},
        .color = {255, 0, 255, 255},
        .radius = 0.5,
        .specular = 0.2,
        .lambert = 0.7,
        .ambient = 0.1,
    };

    Camera *camera     = &scene.camera;
    Vec3    eye_vector = vec3_unit(vec3_sub(camera->orn, camera->pos));
    Vec3    vp_right   = vec3_unit(vec3_cross(eye_vector, (Vec3){0, 1, 0}));
    Vec3    vp_up      = vec3_unit(vec3_cross(vp_right, eye_vector));

    float fov_radians   = deg2rad(camera->fov / 2);
    float aspect_ratio  = WIDTH / HEIGHT;
    float half_width    = tan(fov_radians);
    float half_height   = aspect_ratio * half_width;
    float camera_width  = half_width * 2;
    float camera_height = half_height * 2;
    float pixel_width   = camera_width / (WIDTH - 1);
    float pixel_height  = camera_height / (HEIGHT - 1);

    Color data[HEIGHT * WIDTH];
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            Vec3 xcomp = vec3_scale(vp_right, (x * pixel_width) - half_width);
            Vec3 ycomp = vec3_scale(vp_up, (y * pixel_height) - half_height);
            Ray  ray   = {
                .pos = camera->pos,
                .orn = vec3_unit(vec3_add(vec3_add(eye_vector, xcomp), ycomp)),
            };
            data[WIDTH * y + x] = scene_trace(&scene, ray, 0);
        }
    }

    stbi_write_png("out.png", WIDTH, HEIGHT, 4, data, WIDTH * 4);
    return 0;
}
