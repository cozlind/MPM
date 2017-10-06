using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(Basic))]
[CanEditMultipleObjects]
public class BasicEditor : Editor {
    

    public override void OnInspectorGUI () {
        Basic basic = (Basic)target;

        DrawDefaultInspector ();


        Repaint ();
    }
}
