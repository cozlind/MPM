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
    public float[,] gridMass_w;
    public float[,] gridMass_s;
    public Vector3[,] gridVelocity_w;
    public Vector3[,] gridVelocity_s;
    Matrix4x4 C = Matrix4x4.zero;
    Matrix<double> F_s = Matrix<double>.Build.DenseIdentity (2);
    [Header ("Hydrostatic Axis")]
    public float c0 = 1.1f;
    public float a = -0.5f;
    public float b = 0;
    public float sc = 0.15f;

    public Solver (Body b) {
        body = b;
        gridVelocity_w = new Vector3[22, 22];
        gridVelocity_s = new Vector3[22, 22];
    }

    public void Step(float dt) {
        if (dt == 0.0) return;
        TransferToGrid ();
        UpdateGridMomenta (dt);
    }
    private void TransferToGrid () {
        gridMass_w = new float[22, 22];
        gridMass_s = new float[22, 22];
        //traverse the grid
        for (int i = 0; i < 22; i++) {
            for (int j = 0; j < 22; j++) {
                //water
                gridMass_w[i, j] = 0;
                Vector3 sum_w = Vector3.zero;
                //traverse the particles
                for (int k = 0; k < body.num_w; k++) {
                    Vector3 gridPos = new Vector3 (i * gridSpacing, j * gridSpacing, 0);
                    float w = Kernel (body.pos_w[k] - gridPos);
                    gridMass_w[i, j] += w * body.particleMass_w;
                    sum_w += w * body.particleMass_w * (body.particleVelocity_w[k] + C.MultiplyPoint (gridPos - body.pos_w[k]));
                }
                gridVelocity_w[i, j] = sum_w / gridMass_w[i, j];

                //sand
                gridMass_s[i, j] = 0;
                Vector3 sum_s = Vector3.zero;
                //traverse the particles
                for (int k = 0; k < body.num_s; k++) {
                    Vector3 gridPos = new Vector3 (i * gridSpacing, j * gridSpacing, 0);
                    float w = Kernel (body.pos_s[k] - gridPos);
                    gridMass_s[i, j] += w * body.particleMass_s;
                    sum_s += w * body.particleMass_s * (body.particleVelocity_s[k] + C.MultiplyPoint (gridPos - body.pos_s[k]));
                }
                gridVelocity_s[i, j] = sum_s / gridMass_s[i, j];
            }
        }
    }
    private void UpdateGridMomenta (float dt) {

        for (int i = 0; i < 22; i++) {
            for (int j = 0; j < 22; j++) {
                Vector3 xOld_s = new Vector3 (i * gridSpacing, j * gridSpacing, 0);
                Vector3 xNew_s = xOld_s + dt * gridVelocity_s[i, j];


                //Energy function
                Vector<double> logSingularValues = F_s.Svd (true).S.PointwiseLog ();
                Matrix<double> logStrain = Matrix<double>.Build.DenseOfDiagonalVector(logSingularValues);
                float MU = body.youngsModulus / (1 + body.possionRatio) / 2;
                float LAMBDA = body.youngsModulus * body.possionRatio / (1 + body.possionRatio) / (1 - 2 * body.possionRatio);
                var ElasticPotentialEnergyDensity = MU * logStrain.PointwisePower (2).Trace () + LAMBDA / 2 * logStrain.Trace ();


                //unilateral energy function
                Vector<double> o = DenseVector.OfArray (new double[] { Constants.Sqrt1Over2, Constants.Sqrt1Over2 });//2d
                var u = logSingularValues * o; //var u = logStrain * o;//why should be a scalar??
                var v = (logSingularValues - u * o).L2Norm();
                var f = c0 * v * v * v * v / (1 + Math.Pow (Math.Abs (v), 3));
                var UnilateralEnergyDensity = ElasticPotentialEnergyDensity * H (u, f);


                //weight gradient


                }
        }
    }
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

    private float _CubicBSplineKernel (float x) {
        x = Mathf.Abs (x);
        if (x >= 0 && x < 1) return (x * x * x / 2 - x * x + 2 / 3);
        if (x >= 1 && x < 2) return ((2 - x) * (2 - x) * (2 - x) / 6);
        return 0;
    }
    private float Kernel (Vector3 pos) {
        return _CubicBSplineKernel (pos.x / gridSpacing) * _CubicBSplineKernel (pos.y / gridSpacing) * _CubicBSplineKernel (pos.z / gridSpacing);
    }
    private void UpdateVelocity () {

    }
    private void UpdatePositions () {
        for(int i = 0; i < body.num_w; i++) { 
            body.pos_w[i] += Vector3.one*0.1f;
        }
    }
}
