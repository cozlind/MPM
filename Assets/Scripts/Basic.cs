using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Basic : MonoBehaviour {

    private const float timeStep = 1.0f / 60.0f;
    public float radius = 0.25f;

    public bool drawLines = true;
    public bool drawBoundary = false;

    private List<GameObject> WaterSpheres { get; set; }
    private List<GameObject> SandSpheres { get; set; }
    private List<GameObject> BoundarySpheres { get; set; }

    public Material boundaryMat;
    public Material waterMat;
    public Material sandMat;

    [Header ("Boundary")]
    public float xmin;
    public float xmax;
    public float ymin;
    public float ymax;
    [Header ("Water")]
    public float xmin_w;
    public float xmax_w;
    public float ymin_w;
    public float ymax_w;
    [Header ("Sand")]
    public float xmin_s;
    public float xmax_s;
    public float ymin_s;
    public float ymax_s;

    [Header ("Solver Parameters")]
    public Solver solver;
    public Body body;

    void Start () {
        WaterSpheres = new List<GameObject> ();
        SandSpheres = new List<GameObject> ();
        BoundarySpheres = new List<GameObject> ();

        CreateBoundary (radius);
        CreateWater (radius);
        CreateSand (radius);

        body = new Body (WaterSpheres, SandSpheres, BoundarySpheres);
        solver = new Solver (body);
    }
	
	void FixedUpdate () {
        solver.Step (timeStep);
        UpdateSpheres ();
	}

    private void OnDrawGizmos () {
        if (drawLines) {
            Gizmos.color = Color.black;
            Gizmos.DrawLine (new Vector3 (xmin, ymin, 0), new Vector3 (xmin, ymax, 0));
            Gizmos.DrawLine (new Vector3 (xmin, ymax, 0), new Vector3 (xmax, ymax, 0));
            Gizmos.DrawLine (new Vector3 (xmax, ymin, 0), new Vector3 (xmax, ymax, 0));
            Gizmos.DrawLine (new Vector3 (xmin, ymin, 0), new Vector3 (xmax, ymin, 0));
            Gizmos.color = Color.blue;
            Gizmos.DrawLine (new Vector3 (xmin_w, ymin_w, 0), new Vector3 (xmin_w, ymax_w, 0));
            Gizmos.DrawLine (new Vector3 (xmin_w, ymax_w, 0), new Vector3 (xmax_w, ymax_w, 0));
            Gizmos.DrawLine (new Vector3 (xmax_w, ymin_w, 0), new Vector3 (xmax_w, ymax_w, 0));
            Gizmos.DrawLine (new Vector3 (xmin_w, ymin_w, 0), new Vector3 (xmax_w, ymin_w, 0));
            Gizmos.color = Color.yellow;
            Gizmos.DrawLine (new Vector3 (xmin_s, ymin_s, 0), new Vector3 (xmin_s, ymax_s, 0));
            Gizmos.DrawLine (new Vector3 (xmin_s, ymax_s, 0), new Vector3 (xmax_s, ymax_s, 0));
            Gizmos.DrawLine (new Vector3 (xmax_s, ymin_s, 0), new Vector3 (xmax_s, ymax_s, 0));
            Gizmos.DrawLine (new Vector3 (xmin_s, ymin_s, 0), new Vector3 (xmax_s, ymin_s, 0));
        }
    }

    GameObject CreateSphere(Vector3 pos,float radius,Material mat, bool active=true) {
        GameObject sphere = GameObject.CreatePrimitive (PrimitiveType.Sphere);
        sphere.SetActive (active);
        sphere.transform.parent = transform;
        sphere.transform.position = pos;
        sphere.transform.localScale = Vector3.one * radius;
        sphere.GetComponent<Collider> ().enabled = false;
        sphere.GetComponent<MeshRenderer> ().material = mat;
        return sphere;
    }
    void CreateBoundary (float radius) {
        float half = radius / 2;
        for(float x = xmin - half, y = ymin - half; y <= ymax - half; y+= radius) {
            BoundarySpheres.Add(CreateSphere (new Vector3(x,y,0),radius,boundaryMat, drawBoundary));
        }
        for (float x = xmin - half, y = ymax + half; x <= xmax - half; x+= radius) {
            BoundarySpheres.Add (CreateSphere (new Vector3 (x, y, 0), radius, boundaryMat, drawBoundary));
        }
        for (float x = xmax + half, y = ymax + half; y >= ymin + half; y-=radius) {
            BoundarySpheres.Add (CreateSphere (new Vector3 (x, y, 0), radius, boundaryMat, drawBoundary));
        }
        for (float x = xmax + half, y = ymin - half; x >= xmin + half; x -= radius) {
            BoundarySpheres.Add (CreateSphere ( new Vector3 (x, y, 0), radius, boundaryMat, drawBoundary));
        }
    }
    void CreateWater (float radius) {
        float half = radius / 2;
        for(float y = ymin_w + half; y <= ymax_w - half; y += radius) {
            for(float x = xmin_w + half; x <= xmax_w - half; x += radius) {
                WaterSpheres.Add(CreateSphere ( new Vector3 (x, y, 0), radius, waterMat));
            }
        }
    }
    void CreateSand (float radius) {
        float half = radius / 2;
        for (float y = ymin_s + half; y <= ymax_s - half; y += radius) {
            for (float x = xmin_s + half; x <= xmax_s - half; x += radius) {
                SandSpheres.Add (CreateSphere (new Vector3 (x, y, 0), radius, sandMat));
            }
        }
    }

    public void UpdateSpheres () {

        if (WaterSpheres != null) {
            for (int i = 0; i < WaterSpheres.Count; i++) {
                WaterSpheres[i].transform.position = body.pos_w[i];
            }
        }
        if (SandSpheres != null) {
            for (int i = 0; i < SandSpheres.Count; i++) {
                SandSpheres[i].transform.position = body.pos_s[i];
            }
        }

        if (BoundarySpheres != null) {
            for (int i = 0; i < BoundarySpheres.Count; i++) {
                BoundarySpheres[i].SetActive (drawBoundary);
            }
        }
    }

}
