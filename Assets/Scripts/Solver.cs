using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

[System.Serializable]
public class Solver {


    public float density;
    public Body body { get; set; }
    public float gravity = 9.8f;
    [Header ("Grid")]
    public float gridSpacing = 0.25f;
    public double[,] mass_wg;
    public double[,] mass_sg;
    public Vector<double>[,] Velocity_wg;
    public Vector<double>[,] Velocity_sg;
    Matrix<double> C = Matrix<double>.Build.Dense(2,2,0);
    Matrix<double> F_n_s;
    [Header ("Hydrostatic Axis")]
    public float c0 = 1.1f;
    public float a = -0.5f;
    public float b = 0;
    public float sc = 0.15f;

    public Solver (Body b) {
        body = b;
        Velocity_wg =new Vector<double>[22, 22];
        Velocity_sg = new Vector<double>[22, 22];
        for(int i = 0; i < 22; i++) {
            for(int j = 0; j < 22; j++) {
                Velocity_wg[i, j] = Vector<double>.Build.Dense (2, 0);
                Velocity_sg[i, j] = Vector<double>.Build.Dense (2, 0);
            }
        }
    }

    public void Step(float dt) {
        if (dt == 0.0) return;
        TransferToGrid ();
        UpdateGridMomenta (dt);
    }
    private void TransferToGrid () {
        mass_wg = new double[22, 22];
        mass_sg = new double[22, 22];
        //traverse the grid
        for (int i = 0; i < 22; i++) {
            for (int j = 0; j < 22; j++) {
                //water
                mass_wg[i, j] = 0;
                Vector<double> sum_w = Vector<double>.Build.Dense(2,0);
                //traverse the particles
                for (int k = 0; k < body.num_w; k++) {
                    Vector<double> Pos_g = DenseVector.OfArray (new double[] { i * gridSpacing, j * gridSpacing });
                    Vector<double> Pos_wp = DenseVector.OfArray (new double[] { body.pos_w[k].x, body.pos_w[k].y });
                    Vector<double> Velocity_wp = DenseVector.OfArray (new double[] { body.particleVelocity_w[k].x, body.particleVelocity_w[k].y });

                    double w = Kernel2d (Pos_wp - Pos_g);
                    mass_wg[i, j] += w * body.particleMass_w;
                    sum_w += w * body.particleMass_w * (Velocity_wp + C* (Pos_g - Pos_wp));
                }
                Velocity_wg[i, j] = sum_w / mass_wg[i, j];

                //sand
                mass_sg[i, j] = 0;
                Vector<double> sum_s = Vector<double>.Build.Dense (2, 0);
                //traverse the particles
                for (int k = 0; k < body.num_s; k++) {
                    Vector<double> Pos_g = DenseVector.OfArray (new double[] { i * gridSpacing, j * gridSpacing });
                    Vector<double> Pos_sp = DenseVector.OfArray (new double[] { body.pos_s[k].x, body.pos_s[k].y });
                    Vector<double> Velocity_sp = DenseVector.OfArray (new double[] { body.particleVelocity_s[k].x, body.particleVelocity_s[k].y });

                    double w = Kernel2d (Pos_sp - Pos_g);
                    mass_sg[i, j] += w * body.particleMass_s;
                    sum_s += w * body.particleMass_s * (Velocity_sp + C * (Pos_g - Pos_sp));
                }
                Velocity_sg[i, j] = sum_s / mass_sg[i, j];
            }
        }
    }
    private void UpdateGridMomenta (float dt) {

        for (int i = 0; i < 22; i++) {
            for (int j = 0; j < 22; j++) {
                Vector<double> xOld_s = DenseVector.OfArray(new double[] { i * gridSpacing, j * gridSpacing });
                Vector<double> xNew_s = xOld_s + dt * Velocity_sg[i, j];
                
                //calculate the function paramter F_s
                Matrix<double> sum = Matrix<double>.Build.DenseIdentity (2);
                for (int k = 0; k < body.num_s; k++) {
                    var x_sp = DenseVector.OfArray (new double[] { body.pos_s[k][0], body.pos_s[k][1] });
                    var delW = GradientKernel2d (x_sp - xOld_s);
                    sum += (xNew_s - xOld_s).ToColumnMatrix () * delW.ToRowMatrix ();
                }
                Matrix<double> F_s = sum * F_n_s;




                //Energy function
                var svd = F_s.Svd (true);
                Vector<double> lnSingularValues = svd.S.PointwiseLog10 ();
                Matrix<double> lnStrain = Matrix<double>.Build.DenseOfDiagonalVector(lnSingularValues);
                float MU = body.youngsModulus / (1 + body.possionRatio) / 2;
                float LAMBDA = body.youngsModulus * body.possionRatio / (1 + body.possionRatio) / (1 - 2 * body.possionRatio);
                //var ElasticPotentialEnergyDensity = MU * lnStrain.PointwisePower (2).Trace () + LAMBDA / 2 * lnStrain.Trace ();




                //unilateral energy function
                //Vector<double> o = DenseVector.OfArray (new double[] { Constants.Sqrt1Over2, Constants.Sqrt1Over2 });//2d
                //var u = lnSingularValues * o; //var u = logStrain * o;//why should be a scalar??
                //var v = (lnSingularValues - u * o).L2Norm();
                //var f = c0 * v * v * v * v / (1 + Math.Pow (Math.Abs (v), 3));
                //var UnilateralEnergyDensity = ElasticPotentialEnergyDensity * H (u, f);





                //force and momenta
                var EnergyDensityDerivative = svd.U * (2 * MU * svd.W.Inverse () * lnStrain + LAMBDA * lnStrain.Trace () * svd.W.Inverse ()) * svd.VT;
                double Vp0 = 2 * Constants.Pi * Body.particleRadius * Body.particleRadius;//2d
                var f_s = Vector<double>.Build.Dense (2,  0);
                for (int k = 0; k < body.num_s; k++) {
                    var x_sp = DenseVector.OfArray (new double[] { body.pos_s[k][0], body.pos_s[k][1] });
                    var delW = GradientKernel2d (x_sp - xOld_s);
                    f_s += Vp0 * EnergyDensityDerivative * F_n_s.Transpose () * delW;
                }
                f_s = -f_s;




                //explicit
                //velocity
                f_s += DenseVector.OfArray (new double[] { 0, -1 }) * gravity;
                var v = Velocity_sg[i, j] + f_s * dt / mass_sg[i, j];
                //implicit
                //collision
                //friction
            }
        }
    }
    #region transformation between Math.net and Mathf
    Vector<double> V3ToVec(Vector3 v) {
        return DenseVector.OfArray (new double[] { v[0], v[1], v[2] });
    }
    Vector<double> V2ToVec (Vector3 v) {
        return DenseVector.OfArray (new double[] { v[0], v[1]});
    }
    #endregion

    #region hs() and h()
    private double H(double u,double f) {
        if (u + f < a + sc) return 1;
        if (u + f > b + sc) return 0;
        return _hs ((u + f - a - sc) / (b - a));
    }
    private double _hs (double z) {
        if (z < 0) return 1;
        if (z > 1) return 0;
        return 1 - 10 * Math.Pow (z, 3) + 15 * Math.Pow (z, 4) - 6 * Math.Pow (z, 5);
    }
    #endregion

    #region N() Cubic B Spline Kernel Function
    private double _N (double x) {
        x = Math.Abs (x);
        if (x >= 0 && x < 1) return (x * x * x / 2 - x * x + 2 / 3);
        if (x >= 1 && x < 2) return ((2 - x) * (2 - x) * (2 - x) / 6);
        return 0;
    }
    private double Kernel (Vector<double> pos) {
        pos /= gridSpacing;
        return _N (pos[0]) * _N (pos[1]) * _N (pos[2]);
    }
    private double Kernel2d (Vector<double> pos) {
        pos /= gridSpacing;
        return _N (pos[0]) * _N (pos[1]);
    }
    #endregion

    #region NablaN() Gradient Kernel Function
    private double _DelN (double x) {
        x = Math.Abs (x);
        if (x >= 0 && x < 1) return (1.5f * x * x - 2 * x);
        if (x >= 1 && x < 2) return (-0.5f * x * x + 2 * x - 2);
        return 0;
    }
    private Vector<double> GradientKernel (Vector<double> pos) {
        pos /= gridSpacing;
        var NX = _DelN (pos[0]) * _N (pos[1]) * _N (pos[2]);
        var NY = _N (pos[0]) * _DelN (pos[1]) * _N (pos[2]);
        var NZ = _N (pos[0]) * _N (pos[1]) * _DelN (pos[2]);
        return DenseVector.OfArray(new double[] { NX, NY, NZ });
    }
    private Vector<double> GradientKernel2d (Vector<double> pos) {
        pos /= gridSpacing;
        var NX = _DelN (pos[0]) * _N (pos[1]) * _N (pos[2]);
        var NY = _N (pos[0]) * _DelN (pos[1]) * _N (pos[2]);
        var NZ = _N (pos[0]) * _N (pos[1]) * _DelN (pos[2]);
        return DenseVector.OfArray (new double[] { NX, NY });
    }
    #endregion

    private void UpdateVelocity () {

    }
    private void UpdatePositions () {
        for(int i = 0; i < body.num_w; i++) { 
            body.pos_w[i] += Vector3.one*0.1f;
        }
    }
}
