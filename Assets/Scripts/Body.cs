using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class Body{
    [Header ("Body Info")]
    [SerializeField]
    public static float particleRadius = 0.25f;
    public int num_w;
    public int num_s;
    public int num_b;
    public Vector3[] pos_w;
    public Vector3[] pos_s;
    public Vector3[] pos_b;
    public Vector3[] particleVelocity_w;
    public Vector3[] particleVelocity_s;
    public float particleMass_w =1;
    public float particleMass_s=1;
    [Header ("Physical Info")]
    [Range (0.20f, 0.45f)]
    public float possionRatio = 0.30f;
    [Range (10, 192)]
    public float youngsModulus = 80;

    public Body(List<GameObject> w, List<GameObject> s, List<GameObject> b) {
        num_w = w.Count;
        num_s = s.Count;
        num_b = b.Count;
        pos_w = new Vector3[num_w];
        pos_s = new Vector3[num_s];
        pos_b = new Vector3[num_b];
        particleVelocity_w = new Vector3[num_w];
        particleVelocity_s = new Vector3[num_s];
        for (int i=0;i<num_w;i++) {
            pos_w[i] = w[i].transform.position;
        }
        for (int i = 0; i < num_s; i++) {
            pos_s[i] = s[i].transform.position;
        }
        for (int i = 0; i < num_b; i++) {
            pos_b[i] = b[i].transform.position;
        }
    }

}
