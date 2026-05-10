#include <vector>
#include <functional>
#include <cmath>

#include <include/utils.h>


#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED


// Forest-Ruth 4th order symplectic integrator for N body system
template<int N>
class Forest_Ruth {
    using VectorArray = Array<Vector, N>;

public:
    // constructor
    Forest_Ruth(
        std::function<VectorArray(const VectorArray&, const VectorArray&)> x_dot,
        std::function<VectorArray(const VectorArray&, const VectorArray&)> p_dot,
        VectorArray x0,
        VectorArray p0,
        double t0,
        bool save_results = true)
        :
        x_dot_(x_dot),
        p_dot_(p_dot),
        x_(x0),
        p_(p0),
        t_(t0),
        save_results_(save_results)
    {
        X_.push_back(x_);
        P_.push_back(p_);
        T_.push_back(t_);
    }

    // destructor
    ~Forest_Ruth() = default;

    // get the last time
    double GetCurrentTime() const {
        return t_;
    }

    // get the last position
    VectorArray GetCurrentPosition() const {
        return x_;
    }

    // get the last momentum
    VectorArray GetCurrentMomentum() const {
        return p_;
    }

    // get all times
    std::vector<double> GetTimes() const {
        if (!save_results_) {
            std::cerr << "Warning: GetTimes() called but save_results is false. Returning empty vector." << std::endl;
        }
        return T_;
    }

    // get all computed positions
    std::vector<VectorArray> GetPositions() const {
        if (!save_results_) {
            std::cerr << "Warning: GetPositions() called but save_results is false. Returning empty vector." << std::endl;
        }
        return X_;
    }

    // get all computed momenta
    std::vector<VectorArray> GetMomenta() const {
        if (!save_results_) {
            std::cerr << "Warning: GetMomenta() called but save_results is false. Returning empty vector." << std::endl;
        }
        return P_;
    }

    // make a step of the Forest Ruth algorithm
    void MakeStep(double dt) {
        t_ += dt;
        if (save_results_) {
            T_.push_back(t_);
        }

        std::pair<VectorArray, VectorArray> step = ComputeForestRuthStep(dt, x_, p_, t_);
        x_ = step.first;
        p_ = step.second;
        if (save_results_) {
            X_.push_back(x_);
            P_.push_back(p_);
        }
    }

    // call a loop over MakeStep for K steps at a fixed width
    void MakeSteps(double dt, unsigned int K) {
        for (unsigned int k = 0; k < K; k++) {
            MakeStep(dt);
        }
    }

    // call a loop over MakeStep for a vector of widths
    void MakeStepsVector(const std::vector<double> &dT) {
        for (double dt : dT) {
            MakeStep(dt);
        }
    }


private:
    // function that gives x_dot(x, p)
    std::function<VectorArray(const VectorArray&, const VectorArray&)> x_dot_;
    // function that gives p_dot(x, p)
    std::function<VectorArray(const VectorArray, const VectorArray)> p_dot_;

    // weather to save compuited values or not (set false for long integrations to avoid memory problems)
    bool save_results_;

    // vector that stores the times for which values were calculated
    std::vector<double> T_;
    // vector that stores the computed positions
    std::vector<VectorArray> X_;
    // vector that stores the computed momenta
    std::vector<VectorArray> P_;

    // current time
    double t_;
    // current position
    VectorArray x_;
    // current momentum
    VectorArray p_;

    // algorithm constants
    double c1 = 1/(2-pow(2, 1.0/3));
    double c2 = 1-2*c1;

    // compute a step of the Leapfrog algorithm
    std::pair<VectorArray, VectorArray> ComputeLeapFrogStep(double dt, VectorArray x, VectorArray p, double t) const {
        VectorArray p05 = p + dt/2 * p_dot_(x, p);
        VectorArray x1 = x + dt * x_dot_(x, p05);
        VectorArray p1 = p05 + dt/2 * p_dot_(x1, p05);

        return std::make_pair(x1, p1);
    }

    // compute a step of the Forest-Ruth algorithm
    std::pair<VectorArray, VectorArray> ComputeForestRuthStep(double dt, VectorArray x, VectorArray p, double t) const {

        std::pair<VectorArray, VectorArray> step1 = ComputeLeapFrogStep(c1*dt, x, p, t);
        std::pair<VectorArray, VectorArray> step2 = ComputeLeapFrogStep(c2*dt, step1.first, step1.second, t);
        std::pair<VectorArray, VectorArray> step3 = ComputeLeapFrogStep(c1*dt, step2.first, step2.second, t);

        return step3;
    }

};



#endif //SOLVER_H_INCLUDED

