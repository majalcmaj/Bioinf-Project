"use strict";
var fwd = new Array(), rev = new Array();
fwd[0] = new Array();
fwd[1] = new Array();
rev[0] = new Array();
rev[1] = new Array();
var s1, s2;
var seq1 = "", seq2 = "";

function leftSubMatrix(u_idx1, u_idx2, v_idx1, v_idx2, cols)
// DPA on s1[p1..p2) and s2[q1..q2)
{
    var i, j;

    cols[u_idx1 % 2][v_idx1] = 0.0; // boundary conditions
    for (j = v_idx1 + 1; j <= v_idx2; j++)
        cols[u_idx1 % 2][j] = cols[u_idx1 % 2][j - 1];

    for (i = u_idx1 + 1; i <= u_idx2; i++)     // outer loop
    {
        cols[i % 2][v_idx1] = cols[(i - 1) % 2][v_idx1]; // boundary conditions

        for (j = v_idx1 + 1; j <= v_idx2; j++)  // inner loop
        {
            var diag = cols[(i - 1) % 2][j - 1];
            if (s1.charAt(i - 1) === s2.charAt(j - 1)) diag += 1;

            cols[i % 2][j] = Math.max(diag, cols[(i - 1) % 2][j], cols[i % 2][j - 1])
        }//for j
    }//for i
}//leftSubMatrix

function rightSubMatrix(p1, p2, q1, q2, cols)
// DPA on reverse(s1[p1..p2)) and reverse(s2[q1..q2))
{
    var i, j;

    cols[p2 % 2][q2] = 0.0; // boundary conditions
    for (j = q2 - 1; j >= q1; j--)
        cols[p2 % 2][j] = cols[p2 % 2][j + 1];

    for (i = p2 - 1; i >= p1; i--) {
        cols[i % 2][q2] = cols[(i + 1) % 2][q2];

        for (j = q2 - 1; j >= q1; j--) {
            var diag = cols[(i + 1) % 2][j + 1];
            if (s1.charAt(i) === s2.charAt(j)) diag += 1;

            cols[i % 2][j] = Math.max(diag, cols[(i + 1) % 2][j], cols[i % 2][j + 1]);
        }
    }
}//rightSubMatrix


function _hirschberg(u_idx1, u_idx_2, v_idx_1, v_idx_2)
// _hirschberg s1[p1..p2) with s2[q1..q2)
{
    var mid, i;
    if (u_idx_2 <= u_idx1) // s1 is empty string
        for (i = v_idx_1; i < v_idx_2; i++) {
            seq1 += '-';
            seq2 += s2.charAt(i);
        }

    else if (v_idx_2 <= v_idx_1) // s2 is empty string
        for (i = u_idx1; i < u_idx_2; i++) {
            seq1 += s1.charAt(i);
            seq2 += '-';
        }

    else if (u_idx_2 - u_idx1 === 1) // s1 is one character and s2 is not empty
    {
        var char = s1.charAt(u_idx1), memo = v_idx_1;
        for (i = v_idx_1 + 1; i < v_idx_2; i++) {
            if (s2.charAt(i) === char) {
                memo = i;
            }
        }
        for (i = v_idx_1; i < v_idx_2; i++) {
            if (i === memo) {
                seq1 += char;
            }
            else {
                seq1 += '-';
            }
            seq2 += s2.charAt(i);
        }
    } // a b [l=2] mid=1, a b c [l=3] mid=1, a b c d [l=4] mid=2

    else // p2>p1+1, s1 has at least 2 chars,  divide s1 and conquer
    {
        mid = Math.floor((u_idx1 + u_idx_2) / 2);
        leftSubMatrix(u_idx1, mid, v_idx_1, v_idx_2, fwd);
        rightSubMatrix(mid, u_idx_2, v_idx_1, v_idx_2, rev);
        var v_mid = v_idx_1,
            best = Number.MIN_VALUE;
        for (i = v_idx_1; i <= v_idx_2; i++) // look for cheapest split of s2
        {
            var sum = fwd[mid % 2][i] + rev[mid % 2][i];
            if (sum > best) {
                best = sum;
                v_mid = i;
            }
        }
        _hirschberg(u_idx1, mid, v_idx_1, v_mid);
        _hirschberg(mid, u_idx_2, v_mid, v_idx_2);
    }
}//_hirschberg

s1 = "ACGTACGTACGT";
s2 = "AGTACCTACCGT";
_hirschberg(0, s1.length, 0, s2.length, 1);

console.log(seq1 + '\n' + seq2);