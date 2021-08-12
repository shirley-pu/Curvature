//import {complex, arbitrary, re, im, add, sub, mul, neg, div, pow, sqrt, besselK} from 'math.js';
function getdeltaGdef_old(test){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

//Can use "Document.querySelector()"
    var test=1;
    
    var d0=document.getElementById("d0").value;
    var l=document.getElementByID("l").value;
    var r0=document.getElementByID("r0").value;
    var ka=document.getElementByID("ka").value;
    var kc=document.getElementByID("kc").value;
    var kg=document.getElementByID("kg").value;
    var alpha=document.getElementByID("alpha").value;
    var s=document.getElementByID("s").value;
    var c0=document.getElementByID("c0").value;

    var u0=(d0/l)/2;

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr0 and Kknr0 are 2 x 2 arrays;
    // array[1] is real and imaginary of besselk order 0 
    // array [2] is real [1] and imaginary [2] of besselk order 1 

    function getBessel(x){
        let bessel=[besselK(0,x),besselK(1,x)]
        return bessel
    }

    var Kkpr0=[getBessel(mul(kp,r0))];
    var Kknr0=[getBessel(mul(kn,r0))];
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(kknr0[1],u0)),mul(kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],(-u0))),mul(Kkpr0[0],s)),divi);

    temp=mul(add(mul(Ap,mul(kp2,Kkpr0[0])),mul(An,mul(kn2,Kknr0[0]))),s);
    temp=add(add(temp,mul(Ap,mul(kp,mul(kp2,Kkpr0[1])))),mul(An,mul(kn,mul(kn2,mul(Kknr0[1],u0)))));
    temp=add(temp,re(mul(mul(gamma,u0),s)));
    var res=mul(mul(temp,(-Math.PI)),mul(r0,kc));

    if (abs(imag(res))>imagErrorThreshold) {throw new Exception("Error the free energy has a non zero imaginary part of=" + imag(res));}
    return re(res)/changeParams;
}
