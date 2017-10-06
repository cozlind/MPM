using UnityEngine;
using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

public class TestMathNet : MonoBehaviour
{

    // Use this for initialization
    void Start()
    {
        Matrix<double> A = DenseMatrix.OfArray(new double[,] {
                                                                    {1,1,1,1},
                                                                    {1,2,3,4},
                                                                    {4,3,2,1}});



        var v = DenseVector.OfArray (new double[] { 4, 2 });

        Matrix<double> logStrain = Matrix<double>.Build.DenseOfDiagonalVector (v.PointwiseLog ());
        var o = DenseVector.OfArray (new double[] { Constants.Sqrt1Over2, Constants.Sqrt1Over2 });
        Debug.Log (logStrain);
        Debug.Log (logStrain * o);
        //Debug.Log("Mat: " + A);
    }   

    // Update is called once per frame
    void Update()
    {

    }
}
