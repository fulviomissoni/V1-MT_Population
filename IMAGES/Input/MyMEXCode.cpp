#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#define  PI	3.1415926535897932
#define  PI_2 1.57079632679489661923



using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    ArrayFactory e;
private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
public:
    /* Constructor for the class. */
    MexFunction()
    {
        matlabPtr = getEngine();
    }
    void displayError(std::string errorMessage)
    {
        ArrayFactory factory;
        matlabPtr->feval(u"error", 0, std::vector<Array>({
        factory.createScalar(errorMessage) }));
    }
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        checkArguments(outputs, inputs); //required
        //inputs var
        const double lngt_v  = inputs[0][0];
        const TypedArray<double> v = inputs[1];
        const double lngt_th = inputs[2][0];
        const TypedArray<double> th = inputs[3];
        const double k0      = inputs[4][0];
        const double sx      = inputs[5][0];
        const double sy      = inputs[6][0];
        const double n_frames= inputs[7][0];
        const TypedArray<double> t = inputs[8];
        
        //output var
        TypedArray<double> tmp = e.createArray<double>({sx,sy,n_frames});
        CellArray II = e.createCellArray({ lngt_v,lngt_th });
        for(int f=0;f<lngt_v;f++){
            for(int g=0;g<lngt_th;g++){
                II[f][g] = tmp;
            }
        }
        //nested cycles
        for(int f=0;f<lngt_v;f++){
            for(int g=0;g<lngt_th;g++){
                for(int x=0;x<sx;x++){
                    for(int y=0;y<sy;y++){
                        for(int indt=0;indt<n_frames;indt++){
                            tmp[y][x][indt] = cos(2*PI*k0*((x+1)*cos(th[g]-PI_2) + (y+1)*sin(th[g]-PI_2)) - 2*PI*k0*v[f]*t[indt]);
                        }
                    }
                }
                II[f][g] = tmp;
            }  
        }
        outputs[0] = II;
    }
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        //implementation
        if (inputs.size() != 9) {
            displayError("Nine input required.");
        }
        if (outputs.size() > 2) {
          displayError("Too many outputs specified.");
        }
    }
};