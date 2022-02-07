#include "driver_state.h"
#include <cstring>
#include <limits>
using namespace std;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    long size = width * height;
    state.image_color=new pixel[size];
    state.image_depth=new float[size];
    for(int i=0; i< size; ++i){
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = std::numeric_limits<float>::max();
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.

void render_triangle(driver_state& state){
    int triNum = state.num_vertices/3;
    int v = 0;

    for(int i = 0; i < triNum;++i){
        data_vertex input[3];
        data_geometry output[3];
        for(int j = 0; j < 3; ++j){
            input[j].data = &state.vertex_data[v];
            output[j].data = input[j].data;
            state.vertex_shader(input[j], output[j], state.uniform_data);
            v += state.floats_per_vertex;
        }
        clip_triangle(state, output[0], output[1], output[2], 0);
    }
}

void render(driver_state& state, render_type type)
{
    if(type == render_type::triangle){
        render_triangle(state);
    }
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    float k0[3], k1[3], k2[3];
    float w = state.image_width / 2.0;
    float h = state.image_height / 2.0;
    float k = w;
    float* data = new float[MAX_FLOATS_PER_VERTEX];

    for(int i = 0; i < 3; ++i){
        if(i==1)k=h;
        k0[i] = k * (v0.gl_Position[i] / v0.gl_Position[3] +1) - 0.5;
        k1[i] = k * (v1.gl_Position[i] / v1.gl_Position[3] +1) - 0.5;
        k2[i] = k * (v2.gl_Position[i] / v2.gl_Position[3] +1) - 0.5;
    }

    float x_min = max(min(k0[0],min(k1[0], k2[0])),(float)0);
    float x_max = min(max(k0[0],max(k1[0], k2[0])),(float)state.image_width);
    float y_min = max(min(k0[1],min(k1[1], k2[1])),(float)0);
    float y_max = min(max(k0[1],max(k1[1], k2[1])),(float)state.image_height);
    float area = k0[0]*k1[1]+k1[0]*k2[1]+k2[0]*k0[1]-k0[0]*k2[1]-k1[0]*k0[1]-k2[0]*k1[1];

    for(int j = ceil(y_min); j <y_max; j++){
        for(int i = ceil(x_min); i <x_max; i++){
            float beta = -((k0[1]-k2[1])*i + (k2[0]-k0[0])*j +k0[0]*k2[1] - k2[0]*k0[1])/area;
            float gamma = ((k0[1]-k1[1])*i + (k1[0]-k0[0])*j +k0[0]*k1[1] - k1[0]*k0[1])/area;
            //float alpha = ((k1[1]-k2[1])*i + (k2[0]-k1[0])*j +k1[0]*k2[1] - k2[0]*k1[1])/area;
            float alpha = 1 - beta - gamma;
            int pos = i + j * state.image_width;

            if(alpha<1 && alpha>0 && beta<1 && beta>0 && gamma<1 && gamma>0){
                assert(abs(beta+gamma+alpha-1)<0.01);
                float dep = alpha * k0[2] + beta * k1[2] + gamma * k2[2];
                if(dep < state.image_depth[pos]){
                    data_fragment in{data};
                    data_output out;
                    state.image_depth[pos] = dep;
                    for(int k = 0; k < state.floats_per_vertex; k++){
                        if(state.interp_rules[k]==interp_type::flat){
                            in.data[k] = v0.data[k];
                        }
                        if(state.interp_rules[k]==interp_type::noperspective){
                            in.data[k] = alpha*v0.data[k]+beta*v1.data[k]+gamma*v2.data[k];
                        }
                    }
                    state.fragment_shader(in, out, state.uniform_data);
                    state.image_color[pos]=
                        make_pixel(out.output_color[0]*255, out.output_color[1]*255, out.output_color[2]*255);
                }
            }
        }
    }
    delete[] data;
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

