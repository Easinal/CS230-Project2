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


void calculate_new_vertex(driver_state& state, data_geometry& output, const data_geometry& v1, const data_geometry& v2, int pos, int plane, float* data){
    float smooth_alpha = 0;
    if (pos==1){
		smooth_alpha = (v2.gl_Position[3] - v2.gl_Position[plane]) / (v1.gl_Position[plane] - v1.gl_Position[3] + v2.gl_Position[3] - v2.gl_Position[plane]);
    }
	else{
		smooth_alpha = (-v2.gl_Position[3] - v2.gl_Position[plane]) / (v1.gl_Position[plane] + v1.gl_Position[3] - v2.gl_Position[3] - v2.gl_Position[plane]);
    }
    output.gl_Position = smooth_alpha * v1.gl_Position + (1 - smooth_alpha) * v2.gl_Position;
    float noperspective_alpha = smooth_alpha * v1.gl_Position[3] / (smooth_alpha * v2.gl_Position[3] + (1 - smooth_alpha) * v2.gl_Position[3]);
    for (int i = 0; i < state.floats_per_vertex; i++) {
        if(state.interp_rules[i]==interp_type::flat){
            data[i] = v1.data[i];
        }
        if(state.interp_rules[i]==interp_type::smooth) {
			data[i] = smooth_alpha * v1.data[i] + (1 - smooth_alpha) * v2.data[i];
		}
		if(state.interp_rules[i]==interp_type::noperspective) {
            data[i] = noperspective_alpha * v1.data[i] + (1 - noperspective_alpha) * v2.data[i];
		}
	}
	output.data = data;
}
// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    int flag = 0;
    int plane = 0;
    int pos = 1;
/*
    rasterize_triangle(state, v0, v1, v2);
    return;
*/  
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }else{
        int pos_arr[6]= {1, -1, 1, -1, 1, -1};
        int plane_arr[6]= {2, 2, 1, 1, 0, 0};
        pos = pos_arr[face];
        plane = plane_arr[face];
        
        if(pos*v0.gl_Position[plane]>v0.gl_Position[3])flag+=1;
        if(pos*v1.gl_Position[plane]>v1.gl_Position[3])flag+=2;
        if(pos*v2.gl_Position[plane]>v2.gl_Position[3])flag+=4;
        switch(flag){
            case 0:{
                clip_triangle(state,v0,v1,v2,face+1);
                break;
            }
            case 1:{
                data_geometry output1[3];
                data_geometry output2[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v1.data;
                output1[0].gl_Position = v1.gl_Position;
                output1[1].data = v2.data;
                output1[1].gl_Position = v2.gl_Position;
                calculate_new_vertex(state, output1[2], v1, v0, pos, plane, data1);
                output2[0].data = v2.data;
                output2[0].gl_Position = v2.gl_Position;
                output2[2].data = output1[2].data;
                output2[2].gl_Position = output1[2].gl_Position;
                calculate_new_vertex(state, output2[1], v2, v0, pos, plane, data2);
                
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                clip_triangle(state,output2[0],output2[1],output2[2],face+1);
                break;
            }
            case 2:{
                data_geometry output1[3];
                data_geometry output2[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v0.data;
                output1[0].gl_Position = v0.gl_Position;
                output1[1].data = v2.data;
                output1[1].gl_Position = v2.gl_Position;
                calculate_new_vertex(state, output1[2], v2, v1, pos, plane, data1);
                output2[0].data = v0.data;
                output2[0].gl_Position = v0.gl_Position;
                output2[2].data = output1[2].data;
                output2[2].gl_Position = output1[2].gl_Position;
                calculate_new_vertex(state, output2[1], v0, v1, pos, plane, data2);
                
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                clip_triangle(state,output2[0],output2[1],output2[2],face+1);
                break;
            }
            case 3:{
                data_geometry output1[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v2.data;
                output1[0].gl_Position = v2.gl_Position;
                calculate_new_vertex(state, output1[1], v2, v1, pos, plane, data1);
                calculate_new_vertex(state, output1[2], v2, v0, pos, plane, data2);
                
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                break;
            }
            case 4:{
                data_geometry output1[3];
                data_geometry output2[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v0.data;
                output1[0].gl_Position = v0.gl_Position;
                output1[1].data = v1.data;
                output1[1].gl_Position = v1.gl_Position;
                calculate_new_vertex(state, output1[2], v1, v2, pos, plane, data1);
                output2[0].data = v0.data;
                output2[0].gl_Position = v0.gl_Position;
                output2[2].data = output1[2].data;
                output2[2].gl_Position = output1[2].gl_Position;
                calculate_new_vertex(state, output2[1], v0, v2, pos, plane, data2);
                /*
                cout<<endl<<"face : "<<face<<" case : "<<flag<<endl<< v0.gl_Position<<endl<<v1.gl_Position<<endl<<v2.gl_Position<<endl;
                cout<<endl<<"output1 : "<<endl<< output1[0].gl_Position<<endl<<output1[1].gl_Position<<endl<<output1[2].gl_Position<<endl;
                cout<<endl<<"output2 : "<<endl<< output2[0].gl_Position<<endl<<output2[1].gl_Position<<endl<<output2[2].gl_Position<<endl;
                */
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                clip_triangle(state,output2[0],output2[1],output2[2],face+1);
                break;
            }
            case 5:{
                data_geometry output1[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v1.data;
                output1[0].gl_Position = v1.gl_Position;
                calculate_new_vertex(state, output1[1], v1, v0, pos, plane, data1);
                calculate_new_vertex(state, output1[2], v1, v2, pos, plane, data2);
                
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                break;
            }
            case 6:{
                data_geometry output1[3];
                float* data1 = new float[MAX_FLOATS_PER_VERTEX];
                float* data2 = new float[MAX_FLOATS_PER_VERTEX];

                output1[0].data = v0.data;
                output1[0].gl_Position = v0.gl_Position;
                calculate_new_vertex(state, output1[1], v0, v1, pos, plane, data1);
                calculate_new_vertex(state, output1[2], v0, v2, pos, plane, data2);
                
                clip_triangle(state,output1[0],output1[1],output1[2],face+1);
                break;
            }
            default:
                break;
        }

    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state,v0,v1,v2,face+1);
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
        else k=w;
        k0[i] = k * (v0.gl_Position[i] / v0.gl_Position[3] +1) - 0.5;
        k1[i] = k * (v1.gl_Position[i] / v1.gl_Position[3] +1) - 0.5;
        k2[i] = k * (v2.gl_Position[i] / v2.gl_Position[3] +1) - 0.5;
        //cout<<i<<" : "<<k0[i]<<" "<<k1[i]<<" "<<k2[i]<<endl;
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
                        if(state.interp_rules[k]==interp_type::smooth){							
                            float alpha_p = 0, beta_p = 0, gamma_p = 0, c = 0;

							c = (alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3]);
							alpha_p = alpha / (v0.gl_Position[3] * c);
							beta_p = beta / (v1.gl_Position[3] * c);
							gamma_p = gamma / (v2.gl_Position[3] * c);
                            in.data[k] = alpha_p*v0.data[k]+beta_p*v1.data[k]+gamma_p*v2.data[k];
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

