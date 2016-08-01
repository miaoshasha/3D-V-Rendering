//
//  main.cpp
//  ecs277proj1
//
//  Created by Shasha on 1/16/16.
//  Copyright © 2016 Shasha. All rights reserved.
//

//data infor
//All datasets are binary, storing 8bit voxels for all slices, for all scanlines, and for all voxels.


#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
using std::ostream;
using std::istream;
#include <vector>
using std::vector;


//function declarations
void display();
void menu(int value);
void createmenu();

void display();
void display1();
void display2();
void display3();
void display4();

//More pixel buffers!!
float* PixelBuffer;
float* Buffer1;
float* Buffer2;
float* Buffer3;
float* Buffer4;
float* buffer;
int MainWindow;
int Win1, Win2, Win3, Win4;
int reso = 64;//128; //data size
//data source
//http://www.tc18.org/code_data_set/3D_images.html
int wis = 400;  //window size
int sam = reso*3;  //samples along ray

//int reso = 64; //resolution, 128^3
typedef unsigned char BYTE;

//
void raycast(char c, int flaginter);
//void raycast(int window_size, float (&buffer)[128*128*128], char c, int flaginter)
void read();
void count_datafrequency(float buffer[]);
long getFileSize(FILE *file);
// An unsigned char can store 1 Bytes (8bits) of data (0-255)

class tdp {
public:
    float x, y, z;
    tdp() {x = 0; y = 0; z = 0;}
    tdp(float x0, float y0, float z0) { x = x0; y = y0; z = z0;}
};

inline tdp normvec(tdp p1) {
    double norm = p1.x*p1.x + p1.y*p1.y + p1.z*p1.z;
    norm = sqrt(norm);
    p1.x /= norm;
    p1.y /= norm;
    p1.z /= norm;
    return p1;
}

//transfer function
inline tdp transfunc(float kd) {
    float alstart = kd;
    tdp alpha;
    
    /*/red
    
    if (alstart < 0.1)
        alpha.x = 0.1;//1.5*alstart;
    else if (alstart > 0.7)
        alpha.x = 0.9;//*alstart;
    else
        alpha.x = alstart;//*0.8 + 0.07;


    //blue
    if (alstart < 0.3)
        alpha.z = 0.5*alstart;
    else if (alstart > 0.7)
        alpha.z = 0.1;//7/4*alstart + 0.1-0.3*7/4;
    else
        alpha.z = 0.8*alstart;
    
    
    //green
    alpha.z = alstart;
    
    //isocolor
    */
    if (alstart > 0) {
        alpha.x = alstart;
        alpha.y = 0.5*alstart;
        alpha.z = 0.1*alstart;
    }
    
    
    return alpha;
}

inline float dotpro(tdp v1, tdp v2) {
    float pro;
    pro = v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
    return pro;
}

inline float length(tdp p1) {
    float f;
    f = p1.x*p1.x + p1.y*p1.y + p1.z*p1.z;
    f = sqrtf(f);
    return f;
}


class innpt {
public:
    
    float light[3];//light rgb
    float lov[3];//light source location
    float ka, ks;
    float kd;
    //double Il;
    float K;
    int n;
    float rgbam[3];
    float rgbdis[3];
    
    //double rgbsrc[3];
    
    innpt() {
        light[0] = 1; light[1] = 1; light[2] = 1;//source light rgb
        lov[0] = 300; lov[1] = 300; lov[2] = 300;//light source location
        ka = 0.1;//.01;//ambient
        kd = 0.5;//diffuse
        ks = 0.9;//reflective
        K = 0;//calculated depending on the point
        n = 2;
        tdp I(1,1,1);
        //rgbam[0] = 0; rgbam[1] = 0.2; rgbam[2] = 0.2;   //ambiant light intensity
        //Il = 1;
        rgbdis[0] = 1; rgbdis[1] = 1; rgbdis[2] = 1;    //light source
    }
private:
    void updateK (tdp p) {
        K = sqrt((p.x-lov[0])*(p.x-lov[0])+(p.y-lov[1])*(p.y-lov[1])+(p.z-lov[2])*(p.z-lov[2]));
        K += 20;
    }
    void updatekd(float rda) {
        kd = rda;
    }
    
    
    
public:
    tdp phonc(tdp p, tdp norm, float kd, int vie[3]) {
        updateK(p);
        updatekd(kd);
        tdp rgb;
        //light vector
        tdp lv(lov[0]-p.x, lov[1]-p.y, lov[2]-p.z);
        lv = normvec(lv);
        //view vector
        //---view input
        tdp view(vie[0], vie[1], vie[2]);
        //reflection vector
        tdp reflec;
        if (length(norm) > 0.0001) {
            reflec.x = 2*(norm.x*lv.x)*norm.x-lv.x;
            reflec.y = 2*(norm.y*lv.y)*norm.y-lv.y;
            reflec.z = 2*(norm.z*lv.z)*norm.z-lv.z;
            reflec = normvec(reflec);
            rgb.x = ka + light[0]*(kd*dotpro(lv, norm)+ks*pow(dotpro(view, reflec), n))/K;
            rgb.y = ka + light[1]*(kd*dotpro(lv, norm)+ks*pow(dotpro(view, reflec), n))/K;
            rgb.z = ka + light[2]*(kd*dotpro(lv, norm)+ks*pow(dotpro(view, reflec), n))/K;
            }
            //std::cout << "one point calculated "<< rgb.x << std::endl;
        else {
            rgb.x = 0;
            rgb.y = 0;
            rgb.z = 0;
        }
        
        return rgb;
    }
};



class cube {
    
public:
    int ixsl[2], iysl[2]; //coord of cube, int
    int izsl[2];
    
    //static int num;
    float coord[8*3];
    //this coord is a local coord, depending on the projection
    //float co_cof[8];
    //read from file, scalar for one vertex, eight for a cube
    float kdcolor[8] = {0, 0, 0, 0, 0, 0, 0, 0};   //kd for each vertex
    tdp norm;
    //norm.x = 0; norm.y = 0; norm.z = 0;
    //for cubic interpolation
    float nearpts[6];   //x,y,z, smaller, larger
    float gra[3] = {0, 0, 0};
    //for each vertex
    //kd_d = readkd(x, y, z, c, &buffer);
    float kd = 0; //for the sampling point

    cube (float xray, float yray, float zray){
        ixsl[0] = floor(xray);
        ixsl[1] = ceil(xray);
        iysl[0] = floor(yray);
        iysl[1] = ceil(yray);
        izsl[0] = floor(zray);
        izsl[1] = ceil(zray);
        if (xray == 0)
            ixsl[1] = 1;
        else if (xray == reso -1)
            ixsl[0] = reso - 2;
        if (yray == 0)
            iysl[1] = 1;
        else if (yray == reso -1)
            iysl[0] = reso - 2;
        if (zray == 0)
            izsl[1] = 1;
        else if (zray == reso -1)
            izsl[0] = reso - 2;
    }
    private:
    void coordcrea(int x[2], int y[2], int z[2]){
        /*
         y| /z
         |/
         |---------x
         left hand coord
         */
        coord[0] = x[0];coord[1] = y[0];coord[2] = z[0];  //vertex 1
        coord[3] = x[1];coord[4] = y[0];coord[5] = z[0];
        coord[6] = x[1];coord[7] = y[1];coord[8] = z[0];
        coord[9] = x[0];coord[10] = y[1];coord[11] = z[0];
        
        /*front
         
         4|/------/3
         |       |
         |       | /
         1|-------|/2
         
         */
        
        coord[12] = x[0];coord[13] = y[0];coord[14] = z[1];  //5
        coord[15] = x[1];coord[16] = y[0];coord[17] = z[1];  //6
        coord[18] = x[1];coord[19] = y[1];coord[20] = z[1];  //7
        coord[21] = x[0];coord[22] = y[1];coord[23] = z[1];  //8
        
        
        /*back
         
        8/|-----/|7
         |     / |
         |       |
         5/------/|6
         /      /
         */
    }
    
    
    
    
    float readkd(int x, int y, int z) {
        //coord is in normal grid
        float kd = buffer[reso*reso*z + reso*y + x];
        return kd;
    }
    
    
public:
    void findallkd() {
        // find all the kds for the vertices
        coordcrea(ixsl, iysl, izsl);
        for (int is = 0; is < 8; is++) {
            kdcolor[is] = readkd(coord[is*3+0], coord[is*3+1], coord[is*3+2]);
        }
    }
    void findgrad(int x, int y, int z) {
        //a) middle points   b)edge points
        if (x == reso - 1) {
            nearpts[0] = readkd(x-1, y, z);
            nearpts[1] = readkd(x, y, z);
            gra[0] = (nearpts[1]-nearpts[0]);
        }
        else if (x == 0) {
            nearpts[0] = readkd(x, y, z);
            nearpts[1] = readkd(x+1, y, z);
            gra[0] = (nearpts[1]-nearpts[0]); //x direction
        }
        else {
            nearpts[0] = readkd(x-1, y, z);
            nearpts[1] = readkd(x+1, y, z);
            gra[0] = (nearpts[1]-nearpts[0])/2; //x direction
        }
    
        if (y == reso - 1) {
            nearpts[2] = readkd(x, y-1, z);
            nearpts[3] = readkd(x, y, z);
            gra[1] = (nearpts[3]-nearpts[2]);
        }
        else if (y == 0){
            nearpts[2] = readkd(x, y, z);
            nearpts[3] = readkd(x, y+1, z);
            gra[1] = (nearpts[3]-nearpts[2]); //y direction
        }
        else {
            nearpts[2] = readkd(x, y-1, z);
            nearpts[3] = readkd(x, y+1, z);
            gra[1] = (nearpts[3]-nearpts[2])/2; //y direction
        }
        
        if (z == reso - 1) {
            nearpts[4] = readkd(x, y, z-1);
            nearpts[5] = readkd(x, y, z);
            gra[2] = (nearpts[5]-nearpts[4]);
        }
        else if (z == 0){
            nearpts[4] = readkd(x, y, z);
            nearpts[5] = readkd(x, y, z+1);
            gra[2] = (nearpts[5]-nearpts[4]); //z direction
        }
        else {
            nearpts[4] = readkd(x, y, z-1);
            nearpts[5] = readkd(x, y, z+1);
            gra[2] = (nearpts[5]-nearpts[4])/2; //z direction
        }
        
    }
    
    
    //interpolation
    inline void linCuInterp(float xx[3], int flag) {
        float x = xx[0]-ixsl[0], y = xx[1]-iysl[0], z = xx[2]-izsl[0];
        findallkd();
        //flag 1 linear, 0, cubic
        switch (flag) {
            case 1:
            {
                //linear
                norm.x = -(1-y)*(1-z)*kdcolor[0] + (1-y)*(1-z)*kdcolor[1] - y*(1-z)*kdcolor[2] - (1-y)*z*kdcolor[4] + y*(1-z)*kdcolor[5] + (1-y)*z*kdcolor[5] - y*z*kdcolor[7] + y*z*kdcolor[6];
                
                norm.y = -(1-x)*(1-z)*kdcolor[0] - x*(1-z)*kdcolor[1] + (1-x)*(1-z)*kdcolor[3] - (1-x)*z*kdcolor[4] + x*(1-z)*kdcolor[2] - x*z*kdcolor[5] +(1-x)*z*kdcolor[7] + x*z*kdcolor[6];
                
                norm.z = -(1-x)*(1-y)*kdcolor[0] - x*(1-y)*kdcolor[1] -(1-x)*y*kdcolor[3] + (1-x)*(1-y)*kdcolor[4] -x*y*kdcolor[2] + x*(1-y)*kdcolor[5] + (1-x)*y*kdcolor[7] + x*y*kdcolor[6];
                
                //norm = normvec(norm);
                //1-4, 5-8, ordered as in cube
                
                kd = (1-x)*(1-y)*(1-z)*kdcolor[0] + x*(1-y)*(1-z)*kdcolor[1] + (1-x)*y*(1-z)*kdcolor[3] + (1-x)*(1-y)*z*kdcolor[4] + x*y*(1-z)*kdcolor[2] + x*(1-y)*z*kdcolor[5] + (1-x)*y*z*kdcolor[7] + x*y*z*kdcolor[6];
            }
                break;
            
                
                
            case 0:
            {
                //cubic
                // 1) find the gradient of the points
                //taken care of in findgrad
                findgrad(xx[0], xx[1], xx[2]);
                
                //cubic, bijk's
                float bijk[4*4*4];

                
                // 2) calculate the bigk's
                for (int jz = 0; jz <4; jz++) {
                    for (int jy = 0; jy < 4; jy++) {
                        for (int jx = 0; jx < 4; jx++) {
                            float s = jx/3, t = jy/3, r = jz/3;
                            //trilinear interpolation for the small b's
                            bijk[jz*4*4 + jy*4 + jx] = kdcolor[0] + gra[0]*s + gra[1]*t + gra[2]*r;
                            //(1-s)*(1-t)*(1-r)*kdcolor[0] + s*(1-t)*(1-r)*kdcolor[1] + (1-s)*t*(1-r)*kdcolor[3] + (1-s)*(1-t)*r*kdcolor[4] + s*t*(1-r)*kdcolor[2] + s*(1-t)*r*kdcolor[5] + (1-s)*t*r*kdcolor[7] + s*t*r*kdcolor[6];
                        }
                    }
                }
                
                // 3) normal
                float xs[3] = {x, y, z};
                float bx4[4*3]; //x*4, y*4, z*4
                for (int ijc = 0; ijc < 3; ijc++) {
                    bx4[ijc*4 + 0] = cu0(xs[ijc]);
                    bx4[ijc*4 + 1] = cu1(xs[ijc]);
                    bx4[ijc*4 + 2] = cu2(xs[ijc]);
                    bx4[ijc*4 + 3] = cu3(xs[ijc]);
                }
                float bx3[3*3]; //x*3, y*3, z*3
                for (int ijc = 0; ijc < 3; ijc++) {
                    bx3[ijc*3 + 0] = qu0(xs[ijc]);
                    bx3[ijc*3 + 1] = qu1(xs[ijc]);
                    bx3[ijc*3 + 2] = qu2(xs[ijc]);
                }
                
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            norm.x += (bijk[k*4*4 + j*4 + i+1] - bijk[k*4*4 + j*4 + i])*bx3[i]*bx4[4+j]*bx4[8+k];
                        }
                    }
                }
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 4; k++) {
                            norm.y += (bijk[k*4*4 + (j+1)*4 + i] - bijk[k*4*4 + j*4 + i])*bx4[i]*bx3[3+j]*bx4[8+k];
                        }
                    }
                }
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 3; k++) {
                            norm.z += (bijk[(k+1)*4*4 + j*4 + i] - bijk[k*4*4 + j*4 + i])*bx4[i]*bx4[4+j]*bx3[6+k];
                        }
                    }
                }
                norm.x *= 3;
                norm.y *= 3;
                norm.z *= 3;
                
               // norm = normvec(norm);

                //kd for the sample point
                kd = 0;
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            kd += bijk[k*4*4 + j*4 + i]*bx4[i]*bx4[4+j]*bx4[8+k];
                            //still larger than 1 !!!!
                        }
                    }
                }
                //if (kd > 1)
                    //std::cout<< "Isovalue larger than 1!!! Error!!!" << std::endl;
                
            }//case content end
            break;
                
        }//switch
        
    }//function end
    
    float getkd() {
        return kd;
    }
    
    inline float cu3(float x){
        return x*x*x;
    }
    inline float cu2(float x){
        return 3*x*x*(1-x);
    }
    inline float cu1(float x){
        return 3*(1-x)*(1-x)*x;
    }
    inline float cu0(float x){
        return (1-x)*(1-x)*(1-x);
    }
    inline float qu0(float x) {
        return (1-x)*(1-x);
    }
    inline float qu1(float x) {
        return 2*x*(1-x);
    }
    inline float qu2(float x) {
        return x*x;
    }
    
    
};



BYTE *fileBuf;          // Pointer to our buffered data
int main(int argc, char * argv[]) {
    
    typedef unsigned char BYTE;
    PixelBuffer = new float[wis*2 * wis*2 * 3];
    Buffer1 = new float[wis*wis*3];
    Buffer2 = new float[wis*wis*3];
    Buffer3 = new float[wis*wis*3];
    Buffer4 = new float[wis*wis*3];
    buffer = new float[reso*reso*reso];
    
    //reading the file data into buffer
    read();
    
    //scale to unit window size
    for (int ima = 0; ima <reso*reso*reso; ima++) {
        buffer[ima] /= 255;//maxcor = 250
    }
    //count_datafrequency(buffer);

    char c = 'x';
    int figchoice = 1;
    switch (figchoice) {
            //1. linear; 2. cubic; 3. difference
        case 1: {
            int flag = 1;
            //write rgb to buffer 1~3
            //void raycast(char c, int flaginter)
            raycast('z', flag);
            
            raycast('y', flag);
            
            raycast('x', flag);
            
        }
            break;
        case 2: {
            int flag = 0;
            //void raycast(char c, int flaginter)
            raycast('x', flag);
            
            raycast('y', flag);
            
            raycast('z', flag);
        }
            break;
            
            
        default:
            break;
    }
    int klarge1 = 0, klarge2 = 0, klarge3 = 0;
    for (int i = 0; i < wis*wis; i++) {
        if (Buffer1[i*3] > 0.9)
            klarge1++;
        if (Buffer2[i*3] > 0.9)
            klarge2++;
        if (Buffer3[i*3] > 0.9)
            klarge3++;
    }
    std::cout<< klarge1 <<std::endl;
    std::cout<< klarge2 <<std::endl;
    std::cout<< klarge3 <<std::endl;


    /*
    float m1 = 0;
    float m2 = 0;
    float m3 = 0;
    for (int ima = 0; ima < wis*wis*3; ima++) {
        if (m1 < Buffer1[ima])
            m1 = Buffer1[ima];
        if (m2 < Buffer2[ima])
            m2 = Buffer2[ima];
        if (m3 < Buffer3[ima])
            m3 = Buffer3[ima];
    }
    std::cout << "The max's of each buffer are: " << m1 << std::endl;
    std::cout << m2<< std::endl;
    std::cout << m3<< std::endl;*/
    
    
    //-------------------display--------------------
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    //set window size to 200*200
    glutInitWindowSize(wis*2, wis*2);
    //set window position
    glutInitWindowPosition(200, 200); //it's window position!
    
    //create and set main window title
    MainWindow = glutCreateWindow("Hello Project 1~");
    glClearColor(0, 0, 0, 0); //clears the buffer of OpenGL
    //sets display function
    glutDisplayFunc(display);
    //the 5 parameters are: main window ID, xpos, ypos, width, height
    Win1 = glutCreateSubWindow(MainWindow, 0, 0, wis, wis);
    glutDisplayFunc(display1);
    
    Win2 = glutCreateSubWindow(MainWindow, 0, wis, wis, wis);
    glutDisplayFunc(display2);
    
    Win3 = glutCreateSubWindow(MainWindow, wis, 0, wis, wis);
    glutDisplayFunc(display3);
    
    Win4 = glutCreateSubWindow(MainWindow, wis, wis, wis, wis);
    glutDisplayFunc(display4);

    //create menu here-----------------------------
    createmenu();
    glutMainLoop();//main display loop, will display until terminate
    
    
    
    
    return 0;
}



void display()
{
    //Misc.
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    
    //draws pixel on screen, width and height must match pixel buffer dimension
    glDrawPixels(wis*2, wis*2, GL_RGB, GL_FLOAT, PixelBuffer);
    
    //window refresh
    //glEnd();
    glFlush();
}

//additional display function for subwindows
void display1()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glDrawPixels(wis, wis, GL_RGB, GL_FLOAT, Buffer1);
    glFlush();
}

void display2()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glDrawPixels(wis, wis, GL_RGB, GL_FLOAT, Buffer2);
    glFlush();
}

void display3()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glDrawPixels(wis, wis, GL_RGB, GL_FLOAT, Buffer3);
    glFlush();
}

void display4()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glDrawPixels(wis, wis, GL_RGB, GL_FLOAT, Buffer4);
    glFlush();
}


//make a menu for the operations
void createmenu(){
    
    // Create the menu, this menu becomes the current menu
    int MainMenu = glutCreateMenu(menu);
    glutAddMenuEntry("Change light source location", 1);
    glutAddMenuEntry("Change other phong parameters", 2);
    glutAddMenuEntry("Change object", 3);
    glutAddMenuEntry("Clear scene", 5);
    glutAddMenuEntry("Close window", 0);
    
    // Let the menu respond on the right mouse button
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    
}

void menu(int value) {
    if (value == 1) {
       
    }
    
    
    if (value == 0) {
        glutDestroyWindow(MainWindow);
        exit(0);
    }
    

    //glutPostRedisplay();
    glutSetWindow(Win1);
    glutPostRedisplay();
    glutSetWindow(Win2);
    glutPostRedisplay();
    glutSetWindow(Win3);
    glutPostRedisplay();
    glutSetWindow(Win4);
    glutPostRedisplay();

}




void raycast(char c, int flag) {
    //flaginter: 1 linear, 0 cubic
    float colmap[wis*wis*3];
    
    //float alpha = 0.8;    //opacity
    //count the bright color pixels
    //loop over x,y for the screen coords
    for (int ix = 0; ix < wis; ix++) {
        for (int iy = 0; iy < wis; iy++) {
            tdp scrgb(0, 0, 0);     //color, c-1 is 0.05
            //we have a specific ray, direction determined by direc, e.g. [0 0 1], project to xy plane
            //calculate the RGB for all vertices using the phong model
            //2 sampling points inside every cube
            float x, y, z;
            x = (float)ix/(float)(wis-1)*(reso-1);
            y = (float)iy/(float)(wis-1)*(reso-1);
            
            for (int iz = sam-1; iz > -1; iz--) {
                z = (float)iz/(float)
                (sam-1)*(reso-1);//sam sampling points along the ray
                
                //change to regular grid
                float coord[3];
                switch (c) {
                    case 'x':
                    {
                        coord[0] = y;
                        coord[1] = z;
                        coord[2] = x;
                    }
                        break;
                    case 'y':
                    {
                        coord[0] = x;
                        coord[1] = y;
                        coord[2] = z;
                    }
                        break;
                    case 'z':
                    {
                        coord[0] = z;
                        coord[1] = x;
                        coord[2] = y;
                    }
                        
                }
                
                //find the small cube that contains this point, prepare for interpolation
                cube cube1(coord[0], coord[1], coord[2]);
                //get the kd for all vertices in the cube, prepare for interpolation, done in linCuInterp
                
                //linear or cubic interpolation
                cube1.linCuInterp(coord, flag); //normal and kd calculated and stored in cube
                if (cube1.getkd() > 0) {
                    int meh = 0;
                    //std::cout << mah << std::endl;
                }

                //phong model illumination------------------------
                
                int viewvec[3] = {0, 0, 0};
                if (c == 'x')
                    viewvec[0] = 1;
                else if (c == 'y')
                    viewvec[1] = 1;
                else if (c == 'z')
                    viewvec[2] = 1;
                
                tdp p;
                p.x = coord[0]; p.y = coord[1]; p.z = coord[2];
                innpt innptm1;
                float kd = cube1.getkd();
                
                tdp rgb = innptm1.phonc(p, cube1.norm, kd, viewvec);
                
                tdp alpha = transfunc(kd);
                //parameters, opacity, need a function here
                scrgb.x *= 1-alpha.x;
                scrgb.x += alpha.x*rgb.x;
                scrgb.y *= 1-alpha.y;
                scrgb.y += alpha.y*rgb.y;
                scrgb.z *= 1-alpha.z;
                scrgb.z += alpha.z*rgb.z;
                
                
            }// iz for loop
            
            
            colmap[iy*wis*3 + ix*3 + 0] = scrgb.x;
            colmap[iy*wis*3 + ix*3 + 1] = scrgb.y;
            colmap[iy*wis*3 + ix*3 + 2] = scrgb.z;
            //looped over sample points
            
        }// iy for loop
    }// ix for loop
    //normalize the colors
    float max = 0;
    for (int m = 0; m < wis*wis*3; m++) {
        if (max < colmap[m]) {
            max = colmap[m];
        }
    }
    std::cout << "The max color is "<< max << std::endl;
    
    if (max != 0) {
        for (int m = 0; m < wis*wis*3; m++) {
            colmap[m] /= max;
        }
    }
    
    
    switch (c) {
        case 'x':
        {
            float max = 0;
            for (int ib = 0; ib < wis*wis*3; ib++) {
                Buffer1[ib] = colmap[ib];
            }
            
        }
            break;
        case 'y':
        {
            for (int ib = 0; ib < wis*wis*3; ib++) {
                Buffer2[ib] = colmap[ib];
            }
        }
            break;
        case 'z':
        {
            for (int ib = 0; ib < wis*wis*3; ib++) {
                Buffer3[ib] = colmap[ib];
            }
        }
            break;
    }//switch
    
}


long getFileSize(FILE *file)

{
    long lCurPos, lEndPos;
    
    lCurPos = ftell(file);
    
    fseek(file, 0, 2);
    
    lEndPos = ftell(file);
    
    fseek(file, lCurPos, 0);
    
    return lEndPos;
    
}

void read()
{
    const char *filePath = "/Users/shasha/Documents/courses/ECS277/ecs277proj1/ecs277proj1/fuel.raw";//hydrogenAtom (128), fuel(64), skull (256)
    
    
    FILE *file = NULL;      // File pointer
    
    // Open the file in binary mode using the "rb" format string
    
    // This also checks if the file exists and/or can be opened for reading correctly
    if ((file = fopen(filePath, "rb")) == NULL)
        std::cout << "Could not open specified file" << std::endl;
    else
        std::cout << "File opened successfully" << std::endl;
    
    // Get the size of the file in bytes
    long fileSize = getFileSize(file);
    
    // Allocate space in the buffer for the whole file
    fileBuf = new BYTE[fileSize];
    
    // Read the file in to the buffer
    fread(fileBuf, fileSize, 1, file);
    
      //Now that we have the entire file buffered, we can take a look at some binary infomation
     // Lets take a look in hexadecimal
     
     //for (int i = 0; i < 100; i++)
     
     //printf("%X ", fileBuf[i]);
 
    
    //std::cin.get();
    //std::string fs(fileBuf);
    //float buffer = std::stof(fs);
    for (int i = 0; i < fileSize; i++) {
        buffer[i] = (float) fileBuf[i];
       /* if (i%10000 == 0) {
            std::cout << i << std::endl;
        }*/
    }
    //for (int i = 300000; i < 310000; i++)
        
      //  printf("%f ", buffer[i]);
    //delete[]fileBuf;
    
    fclose(file);   // Almost forgot this
    std::cout << "Reading finished" << std::endl;
    return;
}


void count_datafrequency(float buffer[]) {
    int total = reso*reso*reso;
    int count [total];
    for (int i = 0; i < total; i++) {
        count[i] = 0;
    }
    float value;
    for (int i = 0; i < total; i++) {
        value = buffer[i]*10;
        count[(int) floor(value)]++;
    }
    std::cout<< count << std::endl;
}






