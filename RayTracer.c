/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
 */

#include "utils.h"

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;


struct image* env;
struct image* env1;
struct image* env2;
struct image* env3;
struct image* env4;
struct image* env5;
struct object3D *envion;
struct object3D *envion1;
struct object3D *envion2;
struct object3D *envion3;
struct object3D *envion4;
struct object3D *envion5;


void buildScene(void) {
    // Sets up all objects in the scene. This involves creating each object,
    // defining the transformations needed to shape and position it as
    // desired, specifying the reflectance properties (albedos and colours)
    // and setting up textures where needed.
    // Light sources must be defined, positioned, and their colour defined.
    // All objects must be inserted in the object_list. All light sources
    // must be inserted in the light_list.
    //
    // To create hierarchical objects:
    //   Copy the transform matrix from the parent node to the child, and
    //   apply any required transformations afterwards.
    //
    // NOTE: After setting up the transformations for each object, don't
    //       forget to set up the inverse transform matrix!

    struct object3D *o;
    struct pointLS *l;
    struct point3D p;

    ///////////////////////////////////////
    // TO DO: For Assignment 3 you have to use
    //        the simple scene provided
    //        here, but for Assignment 4 you
    //        *MUST* define your own scene.
    //        Part of your mark will depend
    //        on how nice a scene you
    //        create. Use the simple scene
    //        provided as a sample of how to
    //        define and position objects.
    ///////////////////////////////////////

    // Simple scene for Assignment 3:
    // Insert a couple of objects. A plane and two spheres
    // with some transformations.

    // load environment mapping image
    env=readPPMimage(".\Pics\result_table17.ppm");
    env1=readPPMimage(".\Pics\result_table13.ppm");
    env2=readPPMimage(".\Pics\result_table16.ppm");
    env3=readPPMimage(".\Pics\result_table14.ppm");
    env4=readPPMimage(".\Pics\result_table15.ppm");
    env5=readPPMimage(".\Pics\result_table18.ppm");

    // Let's add a plane
    // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
    //o = newPlane(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);

    //main plane
    o = newPlane(.05, .75, .2, .2, .55, .8, .75, 1, 1, 2);  // Note the plane is highly-reflective (rs=rg=.75) so we
    // should see some reflections if all is done properly.
    // Colour is close to cyan, and currently the plane is
    // completely opaque (alpha=1). The refraction index is
    // meaningless since alpha=1
    Scale(o, 10, 10, 1); // Do a few transforms...
    RotateZ(o, PI / 1.5);
    RotateX(o, PI / 2);
    RotateY(o, -PI /2.25);
    Translate(o, 0, -5, 13);
    invert(&o->T[0][0], &o->Tinv[0][0]); // Very important! compute
    loadTexture(o, ".\Pics\texture_13_by_sirius_sdz-d19qqe1.ppm");
    insertObject(o, &object_list); // Insert into object list


    // Let's add a couple spheres
    o = newSphere(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);
    Scale(o, .8, .8, .8);
    Translate(o, -2, 0, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\images.ppm");
    insertObject(o, &object_list);

    o = newSphere(.05, .95, .15, .15, .75, .95, .55, 1, 1, 6);
    Scale(o, .8, .8, .8);
    Translate(o, 2, 0, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\earth2.ppm");
    insertObject(o, &object_list);

    o = newSphere(.05, .95, .15, .15, .75, .95, .55, 1, 1, 6);
    Scale(o, .8, .8, .8);
    Translate(o, 0, 0, 6);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    //loadTexture(o, ".\Pics\planet_texture_by_kexitt-d7n0it6.ppm");
    loadTexture(o, ".\Pics\planet_texture_3_by_bbbeto-d3btbz3.ppm");
	insertObject(o, &object_list);

    o = newSphere(.05, .95, .15, .15, .75, .95, .55, 1, 1, 6);
    Scale(o, .55, .55, .55);
    Translate(o, 0, 2, 3);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\planet_texture_3_by_bbbeto-d3btbz3.ppm");
    insertObject(o, &object_list);

    o = newSphere(.05, .95, .15, .15, .75, .95, .55, 1, 1, 6);
    Scale(o, .9, .9, .9);
    Translate(o, 3, 2.5, 2);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\Lava Texture (590x700).ppm");
    insertObject(o, &object_list);

    o = newSphere(.05, .95, .55, .55, .75, .95, .55, 1, 1, 6);
    Scale(o, .9, .9, .9);
    Translate(o, 2.7, -.2, 2);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    o->ref = 0;
    o->fra = 1;
    insertObject(o, &object_list);

    o = newCylinder(.05, .95, .75, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .5, .5, 1);
    Translate(o, 0, 1, 0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\7895925634_03d4dfc02e_b.ppm");
    insertObject(o, &object_list);

    o = newCylinder(.05, .95, .75, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .5, .5, 1);
    Translate(o, 0, 0, 0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\7895925634_03d4dfc02e_b.ppm");
    insertObject(o, &object_list);

    o = newCylinder(.05, .95, .75, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .5, .5, 1);
    Translate(o, 0, -1, 0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\7895925634_03d4dfc02e_b.ppm");
    insertObject(o, &object_list);

    o = newCylinder(.05, .95, .75, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .5, .5, 1);
    RotateZ(o, PI/2);
    Translate(o, -1, 0, 0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\7895925634_03d4dfc02e_b.ppm");
    insertObject(o, &object_list);

    o = newCylinder(.05, .95, .75, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .5, .5, 1);
    RotateZ(o, -PI/2);
    Translate(o, 1, 0, 0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\7895925634_03d4dfc02e_b.ppm");
    insertObject(o, &object_list);

    o = newSphere(.05, .95, .95, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .8, .8, .8);
    Translate(o, -3, -.2, 2);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, ".\Pics\36_oilrush_ocean_bottom2_diffusemap.ppm");
    insertObject(o, &object_list);


    o = newSphere(.05, .95, .95, .75, .75, .95, .55, 1, 1, 6);
    Scale(o, .8, .8, .8);
    Translate(o, -2.2, 2, 2);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    o->ref = 1;
    o->fra = 0;
    loadTexture(o, ".\Pics\bf3fc58a-2667-41d3-aa02-58b76ce76486_1000.ppm");
    insertObject(o, &object_list);
    
    // Insert a single point light source.
    p.px = 0;
    p.py = 15.5;
    p.pz = -5.5;
    p.pw = 1;
    l = newPLS(&p, .95, .95, .95);
    
    insertPLS(l, &light_list);

    addAreaLight(4, 4, 0, 0, 1, 0, 15.5, -5.5, 8, 8, .95, .95, .95, &object_list , &light_list);
    addAreaLight(4, 4, 0, 0, 1, 0, 10.5, 2.5, 6, 6, .95, .95, .95, &object_list , &light_list);
    addAreaLight(2, 2, 0, 0, 1, 5, 1.5, -3.5, 4, 5, .95, .95, .95, &object_list , &light_list);    

    // End of simple scene for Assignment 3
    // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
    // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
    //
    // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
    //           the relflectance properties of your objects, and the number and type of light sources
    //           in the scene.
}


void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col) {
    // This function implements the shading model as described in lecture. It takes
    // - A pointer to the first object intersected by the ray (to get the colour properties)
    // - The coordinates of the intersection point (in world coordinates)
    // - The normal at the point
    // - The ray (needed to determine the reflection direction to use for the global component, as well as for
    //   the Phong specular component)
    // - The current racursion depth
    // - The (a,b) texture coordinates (meaningless unless texture is enabled)
    //
    // Returns:
    // - The colour for this ray (using the col pointer)
    //

    struct colourRGB tmp_col; // Accumulator for colour components
    double R, G, B; // Colour for the object in R G and B

    // This will hold the colour as we process all the components of
    // the Phong illumination model
    tmp_col.R = 0;
    tmp_col.G = 0;
    tmp_col.B = 0;

    if (obj->texImg == NULL) // Not textured, use object colour
    {
        R = obj->col.R*obj->alpha;
        G = obj->col.G*obj->alpha;
        B = obj->col.B*obj->alpha;
    } else {
        // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
        // for the object. Note that we will use textures also for Photon Mapping.
        obj->textureMap(obj->texImg, a, b, &R, &G, &B);
    }

    //////////////////////////////////////////////////////////////
    // TO DO: Implement this function. Refer to the notes for
    // details about the shading model.
    //////////////////////////////////////////////////////////////

    /************************************************************/
    pointLS *ligsrc = light_list;


    int num = 0;
/*--------------------- this is the dot light ---------------------*/

    while (ligsrc){

      //calculate s's information
      point3D s = ligsrc->p0;
      subVectors(p, &s);
      normalize(&s);

      //calculate r's information 
      point3D r = *n;
      double n_dot_s = dot(n, &s);
      r.px = n_dot_s * r.px;
      r.py = n_dot_s * r.py;
      r.pz = n_dot_s * r.pz;
      subVectors(&s, &r);
      normalize(&r);

      //calculate b's information
      //b points to direction of camera
      //b is the inverse way of ray 
      point3D b;
      b.px = -ray->d.px;
      b.py = -ray->d.py;
      b.pz = -ray->d.pz;
      b.pw = 1;
      normalize(&b);

      //calculate diffuse light component
      double diffuse = max(0, n_dot_s);
      
      //calculate specular light component
      double r_dot_b = dot(&r, &b);
      double specular = pow (max(0, r_dot_b),obj->shinyness);

      //setup the ray from light to intersection
      ray3D light;
      light.d = s;
      light.p0 = *p;

      //find the first hit object of the light ray
      double templ, tempa, tempb;
      object3D *firsthit;
      point3D tempp, tempn;
      findFirstHit(&light, &templ, obj, &firsthit, &tempp, &tempn, &tempa, &tempb);

      //if the light ray hit nothing or hit a light source
      if (templ < 0 || firsthit->isLightSource){
        if (obj->fra){
            tmp_col.R += obj->alb.rs * ligsrc->col.R * specular; 
            tmp_col.G += obj->alb.rs * ligsrc->col.G * specular;
            tmp_col.B += obj->alb.rs * ligsrc->col.B * specular;
        }else{
            tmp_col.R += obj->alb.rd * R * diffuse + obj->alb.rs * ligsrc->col.R * specular; 
            tmp_col.G += obj->alb.rd * G * diffuse + obj->alb.rs * ligsrc->col.G * specular;
            tmp_col.B += obj->alb.rd * B * diffuse + obj->alb.rs * ligsrc->col.B * specular;
        }
      }

      num++;
      ligsrc = ligsrc->next;

    }

    if(obj->fra){
      col->R += tmp_col.R * (1/(depth + 1)) / num;
      col->G += tmp_col.G * (1/(depth + 1)) / num;
      col->B += tmp_col.B * (1/(depth + 1)) / num;
    }else{

    // for ambient colour
    tmp_col.R = R * obj->alb.ra + tmp_col.R/num;
    tmp_col.G = G * obj->alb.ra + tmp_col.G/num;
    tmp_col.B = B * obj->alb.ra + tmp_col.B/num;

/*---------------------above is the dot light ---------------------*/


    // decreasing specular effect
    col->R =  tmp_col.R * (1/(depth + 1));
    col->G =  tmp_col.G * (1/(depth + 1));
    col->B =  tmp_col.B * (1/(depth + 1));
}
    //update col
    /************************************************************/

    // Be sure to update 'col' with the final colour computed here!
    return;

}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b) {
    // Find the closest intersection between the ray and any objects in the scene.
    // It returns:
    //   - The lambda at the intersection (or < 0 if no intersection)
    //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
    //   - The location of the intersection point (in p)
    //   - The normal at the intersection point (in n)
    //
    // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
    // return a self-intersection due to numerical errors for recursive raytrace calls.
    //

    /////////////////////////////////////////////////////////////
    // TO DO: Implement this function. See the notes for
    // reference of what to do in here
    /////////////////////////////////////////////////////////////

    /**********************************************************/
    double tempv = -1.0;
    struct object3D * ObjPtr = object_list;
    //if not NULL go to the next
    while (ObjPtr) {

        point3D tempp, tempn;
        double templ;
        double tempa = 0;
        double tempb = 0;

        //calculate intersection to every object in the list
        if (ObjPtr->intersect(ObjPtr, ray, &templ, &tempp, &tempn, &tempa, &tempb)) {
            //update at first time or meet a closer object
            if ((tempv < 0 || templ < tempv) && ObjPtr != Os) {
                tempv = templ;
                *obj = ObjPtr;
                *p = tempp;
                *n = tempn;
                *a = tempa;
                *b = tempb;
            }
        }
        ObjPtr = ObjPtr->next;
    }

    //update lamda
    if (tempv <= 0) {
        *lambda = -1.0;
    } else {    
        *lambda = tempv;
    }

    return;
    /**********************************************************/
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os) {
    // Ray-Tracing function. It finds the closest intersection between
    // the ray and any scene objects, calls the shading function to
    // determine the colour at this intersection, and returns the
    // colour.
    //
    // Os is needed for recursive calls to ensure that findFirstHit will
    // not simply return a self-intersection due to numerical
    // errors. For the top level call, Os should be NULL. And thereafter
    // it will correspond to the object from which the recursive
    // ray originates.
    //

    double lambda; // Lambda at intersection
    double a, b; // Texture coordinates
    struct object3D *obj; // Pointer to object at intersection
    struct point3D p; // Intersection point
    struct point3D n; // Normal at intersection
    struct colourRGB I = {0, 0, 0}; // Colour returned by shading function

    if (depth > MAX_DEPTH) // Max recursion depth reached. Return invalid colour.
    {
        col->R = -1;
        col->G = -1;
        col->B = -1;
        return;
    }

    ///////////////////////////////////////////////////////
    // TO DO: Complete this function. Refer to the notes
    // if you are unsure what to do here.
    ///////////////////////////////////////////////////////

    /****************************************************/
    colourRGB Phong = {0, 0, 0}; // Phong illumination
    colourRGB Global = {0, 0, 0}; // Global specular illumination

    if ((Os!= NULL) && (Os->fra)){
        findFirstHit(ray, &lambda, NULL, &obj, &p, &n, &a, &b);
    }else{
        findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
    }
    
    //valid color && update color by phong and global speculation
    if (lambda > 0) {

        //get the phong illum for current intersection
        rtShade(obj, &p, &n, ray, depth, a, b, &Phong);
        I.R += Phong.R;
        I.G += Phong.G;
        I.B += Phong.B;


/********************* mirror reflection ************************/

        ray3D r;

        double d_dot_n = dot(&ray->d, &n);
        point3D dnn;
        dnn.px = n.px * 2 * d_dot_n;
        dnn.py = n.py * 2 * d_dot_n;
        dnn.pz = n.pz * 2 * d_dot_n;

        r.d.px = ray->d.px;
        r.d.py = ray->d.py;
        r.d.pz = ray->d.pz;
        r.d.pw = ray->d.pw;

        subVectors(&dnn, &r.d);
        normalize(&r.d);
        
        r.p0.px = p.px;
        r.p0.py = p.py;
        r.p0.pz = p.pz;
        r.p0.pw = p.pw;

        if (obj->ref){
/********************** glossy reflection ************************/
          int idx;
          for (idx = 0; idx < 10; idx++){
          // super sampling to denoise
            point3D * u = cross(&r.d, &n);
            point3D * v = cross(&r.d, u);
            colourRGB tempGlobal = {0,0,0};
            double theta = 2 * PI  * (drand48() * 0.3);
            double phi = 2 * PI * (drand48() * 0.3);
            double dx = sin(theta) * cos(phi);
            double dy = sin(theta) * sin(phi);
            double dz = cos(theta);

            r.d.px = dx * u->px + dy * v->px + dz * r.d.px;
            r.d.py = dx * u->py + dy * v->py + dz * r.d.py;
            r.d.pz = dx * u->pz + dy * v->pz + dz * r.d.pz;
            normalize(&r.d); 

            rayTrace(&r, depth + 1, &tempGlobal, obj);


            if (tempGlobal.R > 0 && tempGlobal.G >0 && tempGlobal.B > 0) {
                Global.R += tempGlobal.R;
                Global.G += tempGlobal.G;
                Global.B += tempGlobal.B;
            }
        }

        Global.R /=10;
        Global.G /=10;
        Global.B /=10;
/********************** glossy reflection ************************/
        }
        else{
/********************* regular reflection ************************/          
          rayTrace(&r, depth + 1, &Global, obj);
/********************* regular reflection ************************/
        }

        if(obj->fra){

          colourRGB tempGlobal = {0, 0, 0};
          ray3D refrac;
          double costheta = dot(&(r.d),&n);
          double sintheta = sqrt(1 - pow(costheta,2));
          double theta = acos(costheta);
          double phi = theta*3/4;
          double sinphi = sin(phi);
          double cosphi = cos(phi);

          //add a portion of nornal to avoiding intersect itself
          refrac.p0.px = r.p0.px - n.px*0.01;
          refrac.p0.py = r.p0.py - n.py*0.01;
          refrac.p0.pz = r.p0.pz - n.pz*0.01;
          refrac.p0.pw = r.p0.pw;
            
          //by refraction equation
          refrac.d.px = (sinphi / sintheta) * (ray->d.px + n.px * costheta) - n.px * cosphi;
          refrac.d.py = (sinphi / sintheta) * (ray->d.py + n.py * costheta) - n.py * cosphi;
          refrac.d.pz = (sinphi / sintheta) * (ray->d.pz + n.pz * costheta) - n.pz * cosphi;
          refrac.d.pw = ray->d.pw; 
            
          rayTrace(&refrac, depth + 1, &tempGlobal, obj);

          if (tempGlobal.R > 0 && tempGlobal.G >0 && tempGlobal.B > 0) {
              Global.R += tempGlobal.R;
              Global.G += tempGlobal.G;
              Global.B += tempGlobal.B;
          }

        }

        if (Global.R > 0 && Global.G >0 && Global.B > 0) {
            I.R += obj->alb.rg * Global.R;
            I.G += obj->alb.rg * Global.G;
            I.B += obj->alb.rg * Global.B;
        }

        //erase overflow colors
        if (I.R > 1) {
          I.R = 1;
        }
        if (I.G > 1) {
          I.G = 1;
        }
        if (I.B > 1) {
          I.B = 1;
        }

        col->R = I.R;
        col->G = I.G;
        col->B = I.B;
        return;

    } else { //invalid color

        int index = -1;
        double a, b, R, G, B;

        // compute a,b vector
        cubeMap(&index, &a, &b, ray->d.px, ray->d.py, ray->d.pz);

         if(index == 0)
           texMap(env, a, b, &R, &G, &B); 
         if(index == 1)
          texMap(env1, a, b, &R, &G, &B); 
        if(index == 2)
          texMap(env2, a, b, &R, &G, &B);
        if(index == 3)
          texMap(env3, a, b, &R, &G, &B);
        if(index == 4)
          texMap(env4, a, b, &R, &G, &B);
        if(index == 5)
          texMap(env5, a, b, &R, &G, &B);

        col->R = R;
        col->G = G;
        col->B = B;

        return;
    }

    /****************************************************/
}

int main(int argc, char *argv[]) {
    // Main function for the raytracer. Parses input parameters,
    // sets up the initial blank image, and calls the functions
    // that set up the scene and do the raytracing.
    struct image *im; // Will hold the raytraced image
    struct view *cam; // Camera and view for this scene
    int sx; // Size of the raytraced image
    int antialiasing; // Flag to determine whether antialiaing is enabled or disabled
    char output_name[1024]; // Name of the output file for the raytraced .ppm image
    struct point3D e; // Camera view parameters 'e', 'g', and 'up'
    struct point3D g;
    struct point3D up;
    double du, dv; // Increase along u and v directions for pixel coordinates
    struct point3D pc, d; // Point structures to keep the coordinates of a pixel and
    // the direction or a ray
    struct ray3D *ray; // Structure to keep the ray from e to a pixel
    struct colourRGB col; // Return colour for raytraced pixels
    struct colourRGB background; // Background colour
    int i, j; // Counters for pixel coordinates
    unsigned char *rgbIm;

    if (argc < 5) {
        fprintf(stderr, "RayTracer: Can not parse input parameters\n");
        fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
        fprintf(stderr, "   size = Image size (both along x and y)\n");
        fprintf(stderr, "   rec_depth = Recursion depth\n");
        fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
        fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }
    sx = atoi(argv[1]);
    MAX_DEPTH = atoi(argv[2]);
    if (atoi(argv[3]) == 0) antialiasing = 0;
    else antialiasing = 1;
    strcpy(&output_name[0], argv[4]);

    fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
    fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
    if (!antialiasing) fprintf(stderr, "Antialising is off\n");
    else fprintf(stderr, "Antialising is on\n");
    fprintf(stderr, "Output file name: %s\n", output_name);

    object_list = NULL;
    light_list = NULL;

    // Allocate memory for the new image
    im = newImage(sx, sx);
    if (!im) {
        fprintf(stderr, "Unable to allocate memory for raytraced image\n");
        exit(0);
    } else rgbIm = (unsigned char *) im->rgbdata;

    ///////////////////////////////////////////////////
    // TO DO: You will need to implement several of the
    //        functions below. For Assignment 3, you can use
    //        the simple scene already provided. But
    //        for Assignment 4 you need to create your own
    //        *interesting* scene.
    ///////////////////////////////////////////////////
    buildScene(); // Create a scene. This defines all the
    // objects in the world of the raytracer

    //////////////////////////////////////////
    // TO DO: For Assignment 3 you can use the setup
    //        already provided here. For Assignment 4
    //        you may want to move the camera
    //        and change the view parameters
    //        to suit your scene.
    //////////////////////////////////////////

    // Mind the homogeneous coordinate w of all vectors below. DO NOT
    // forget to set it to 1, or you'll get junk out of the
    // geometric transformations later on.

    // Camera center is at (0,0,-1)
    e.px = 0;
    e.py = 0;
    e.pz = -3;
    e.pw = 1;

    // To define the gaze vector, we choose a point 'pc' in the scene that
    // the camera is looking at, and do the vector subtraction pc-e.
    // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
    g.px = 0;
    g.py = 0;
    g.pz = 1;
    g.pw = 1;

    // Define the 'up' vector to be the Y axis
    up.px = 0;
    up.py = 1;
    up.pz = 0;
    up.pw = 1;

    // Set up view with given the above vectors, a 4x4 window,
    // and a focal length of -1 (why? where is the image plane?)
    // Note that the top-left corner of the window is at (-2, 2)
    // in camera coordinates.
    cam = setupView(&e, &g, &up, -3, -2, 2, 4);

    if (cam == NULL) {
        fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
        cleanup(object_list, light_list);
        deleteImage(im);
        exit(0);
    }

    // Set up background colour here
    background.R = 0;
    background.G = 0;
    background.B = 0;

    // Do the raytracing
    //////////////////////////////////////////////////////
    // TO DO: You will need code here to do the raytracing
    //        for each pixel in the image. Refer to the
    //        lecture notes, in particular, to the
    //        raytracing pseudocode, for details on what
    //        to do here. Make sure you undersand the
    //        overall procedure of raytracing for a single
    //        pixel.
    //////////////////////////////////////////////////////
    du = cam->wsize / (sx - 1); // du and dv. In the notes in terms of wl and wr, wt and wb,
    dv = -cam->wsize / (sx - 1); // here we use wl, wt, and wsize. du=dv since the image is
    // and dv is negative since y increases downward in pixel
    // coordinates and upward in camera coordinates.

    fprintf(stderr, "View parameters:\n");
    fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
    fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
    printmatrix(cam->C2W);
    fprintf(stderr, "World to camera conversion matrix\n");
    printmatrix(cam->W2C);
    fprintf(stderr, "\n");

    fprintf(stderr, "Rendering row: ");

    #pragma omp parallel for
    for (j = 0; j < sx; j++) // For each of the pixels in the image
    {
        fprintf(stderr, "%d/%d, ", j, sx);

        #pragma omp parallel for
        for (i = 0; i < sx; i++) {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing!
            ///////////////////////////////////////////////////////////////////

            /****************************************************************/

            /* -----------------this is simple ray tracing------------------*/ 
            //camera origin 
            // point3D p0 = {0, 0, 0, 1};

            // //camera direction 
            // point3D d = {cam->wl + (i + 0.5) * du, cam->wt + (j + 0.5) * dv, -3, 0};
            
            // //color of the graph
            // colourRGB col = {0, 0, 0};

            // //create ray from camera position
            // ray3D ray;
            // matVecMult(cam->C2W, &d);
            // matVecMult(cam->C2W, &p0);
            // normalize(&d);
            // ray.p0 = p0;
            // ray.d = d;

            // rayTrace(&ray, 0, &col, NULL);
            
            /* ----------------above is simple ray tracing------------------*/
            
            /*-----------------this is anti-aliasing tracing----------------*/
            colourRGB col = {0, 0, 0};

            int xway = 0;
            int yway = 0;
            for (xway = 0; xway <3; xway++) {
              for (yway = 0; yway <3; yway++) {

                point3D p0 = {0, 0, 0, 1};

                double xran = drand48()/3.0;
                double yran = drand48()/3.0;

                double xkey = xran + xway/3.0;
                double ykey = yran + yway/3.0;
    
                point3D d = {cam->wl + (i + xkey) * du, cam->wt + (j + ykey) * dv, -3, 0};
    
                //create ray from camera position
                ray3D ray;
                matVecMult(cam->C2W, &d);
                matVecMult(cam->C2W, &p0);
                normalize(&d);
                ray.p0 = p0;
                ray.d = d;

                colourRGB tempc;
                rayTrace(&ray, 0, &tempc, NULL);

                if (tempc.R > 0 || tempc.G>0 || tempc.B >0){
                  col.R += tempc.R;
                  col.G += tempc.G;
                  col.B += tempc.B;
                }

              }
            }
            

            col.R /= 9.0;
            col.G /= 9.0;
            col.B /= 9.0;

          /*------------------above is anti-aliasing tracing----------------*/


            //update pixels color
            if (col.R >= 0 || col.G >= 0 || col.B >= 0) { //color of object
                rgbIm[(j * sx + i)* 3 + 0] = (unsigned char)(col.R * 255.0);
                rgbIm[(j * sx + i)* 3 + 1] = (unsigned char)(col.G * 255.0);
                rgbIm[(j * sx + i)* 3 + 2] = (unsigned char)(col.B * 255.0);
            } else { //color of background
                rgbIm[(j * sx + i)* 3 + 0] = (unsigned char)(background.R * 255.0);
                rgbIm[(j * sx + i)* 3 + 1] = (unsigned char)(background.G * 255.0);
                rgbIm[(j * sx + i)* 3 + 2] = (unsigned char)(background.B * 255.0);
            }

        } // end for i
    } // end for j

    
    fprintf(stderr, "\nDone!\n");

    // Output rendered image
    imageOutput(im, output_name);

    // Exit section. Clean up and return.
    cleanup(object_list, light_list); // Object and light lists
    deleteImage(im); // Rendered image
    free(cam); // camera view
    exit(0);
}


        
