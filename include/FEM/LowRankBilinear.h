#ifndef model_LowRankBilinear_H_
#define model_LowRankBilinear_H_
#include <Eigen/Dense>
namespace model
{
    // Optional positive-semidefinite contribution F^T F, separate from sparse triplets.
    template<class Form>
    Eigen::MatrixXd lowRankBilinearRows(const Form& form)
    {
        if constexpr(requires { form.lowRankPenaltyRows(); })
        {
            return form.lowRankPenaltyRows();
        }
        else
        {
            return Eigen::MatrixXd(0,form.field().gSize());
        }
    }
}
#endif
