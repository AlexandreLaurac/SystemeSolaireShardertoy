//////////////////////// Constantes ////////////////////////

#define R iResolution.xy  // résolution de l'écran
#define PI 3.14159265359  // constante PI

// Caméra
#define CP 50.       // distance de la caméra à l'origine en coordonnées sphériques
#define PHI_DEG 70.  // angle de la caméra avec l'axe des y en coordonnées sphériques          

// Algorithme de ray-marching
#define DIST_MIN 0.001          // distance minimale aux objets
#define DIST_MAX 500.0          // distance maximale aux objets
#define RAY_MARCH_STEPS 1000    // nombre maximal de pas pour l'algorithme
#define RAY_MARCH_PRECI 0.0001  // objet considéré comme atteint si l'écart à la surface est inférieure à cette valeur

// Planètes
#define NB_ASTRES 10      // nombre d'astres de la scène
#define RAYON_SOLEIL 1.0  // rayon du soleil - valeur plutôt arbitraire (voir rapport)
#define FT 0.05           // fréquence de révolution de la Terre autour du Soleil (365.256 jours)
#define RT 1.             // Rayon Terrestre (6378 km)
#define UA 10.            // Unité Astronomique (distance Terre - Soleil : 149 600 000 km)

// Illumination
#define COEF_AMB 0.15  // coefficient pondérateur du terme ambiant du modèle de Phong - fixe la couleur des zones d'ombre


//////////////////////// Structures de données ////////////////////////

// Structure rayon
struct Ray
{
    vec3 o ; // origine
    vec3 d ; // direction
} ;

// Structure de données géométriques d'un astre
struct DonneesGeoAstre
{
    float rayAstre ; // rayon de l'astre
    float rayTraj ;  // rayon de la trajectoire de l'astre
    float freqRev ;  // fréquence de révolution de l'astre (autour d'un autre astre)
} ;

// Structure Astre
struct Astre
{
    vec3 position ;  // vecteur position du centre de l'astre
    float dist ;     // dist signée d'un point quelconque à l'astre
    vec3 couleur ;   // couleur initiale de l'astre
    bool type ;      // type d'astre ; si true : soleil, si false : tout autre astre
} ;


/////////////////// Fonctions de distance signées ///////////////////

// Fonction de distance signée à une sphère - Inigo Quilez
float sdSphere (vec3 p, float s)
{
    return length(p)-s ;
}

// Fonction de distance signée à un cylindre "vertical" - Inigo Quilez
float sdCappedCylinder (vec3 p, float h, float r)
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h) ;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0)) ;
}

// Fonction de distance signée à un cylindre fourni par deux vecteurs a et b - Inigo Quilez
float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
    vec3  ba = b - a;
    vec3  pa = p - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);
    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    
    return sign(d)*sqrt(abs(d))/baba;
}


/////////////////// Fonctions d'opérations ensemblistes ///////////////////

// Fonction d'union de deux astres
Astre opUnion (Astre A1, Astre A2)
{
    if (A1.dist < A2.dist) {
        return A1 ;
    }
    else {
        return A2 ;
    }
}

// Fonction de différence de deux astres (astre A1 soustrait à l'astre A2)
Astre opDifference (Astre A1, Astre A2)
{
    if (-A1.dist > A2.dist) {
        return Astre(A2.position, -A1.dist, vec3(0.), false) ;
    }
    else {
        return A2 ;
    }
}

// Fonction de différence de deux objets bis (objet A1 fourni par sa distance signée soustrait à l'astre A2)
Astre opDifference (float distSig1, Astre A2)
{
    if (-distSig1 > A2.dist) {
        return Astre(A2.position, -distSig1, A2.couleur, A2.type) ;
    }
    else {
        return A2 ;
    }
}

// Fonction retournant l' "astre minimal" à partir du tableau des astres
Astre calculAstreMin (Astre a[NB_ASTRES])
{
    // Initialisation
    int i = 0 ;
    Astre astreMin = a[0] ;
    
    // Iteration
    for (i=1 ; i<NB_ASTRES ; i++) {
        if (a[i].dist < astreMin.dist) {
            astreMin = a[i] ; } }

    // Valeur de retour
    return astreMin ;
}


/////////////////// Fonctions de transformations de l'espace ///////////////////

// Matrice de rotation autour de l'axe Y - Nathan Vaughn, dit "inspirnathan"
mat3 rotateY (float theta)
{
    float c = cos(theta) ;
    float s = sin(theta) ;
    return mat3(
            vec3(c, 0, s),
            vec3(0, 1, 0),
            vec3(-s, 0, c)
        ) ;
}


///////////// Fonctions pseudo-aléatoires (code fourni en cours) /////////////

float hash11(float p)
{
    return fract(sin(p*727.1)*43758.5453123) ;
}

float hash12(vec2 p)
{
    return fract (sin(dot(p,vec2(12.9898,78.233)))*43758.5453123) ;
}

float hash13(vec3 p)
{
    return fract(sin(dot(p,vec3(12.9898,78.233,45.5432)))*43758.5453123) ;
}

float vnoise(in vec3 x)
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);

    return mix(
            mix(mix(hash13(p+vec3(0.,0.,0.)),hash13(p+vec3(1.,0.,0.)),f.x),
                  mix(hash13(p+vec3(0.,1.,0.)),hash13(p+vec3(1.,1.,0.)),f.x),f.y),
            mix(mix(hash13(p+vec3(0.,0.,1.)),hash13(p+vec3(1.,0.,1.)),f.x),
                  mix(hash13(p+vec3(0.,1.,1.)),hash13(p+vec3(1.,1.,1.)),f.x),f.y),f.z);
}

float fnoise(in vec3 p, in float amplitude, in float frequency, in float persistence, in int nboctaves)
{
    float a = amplitude ;
    float f = frequency ;
    float n = 0.0 ;

    for (int i=0 ; i<nboctaves ; i++)
    {
        n = n + a*vnoise(p*f) ;
        f = 2. * f ;
        a = a * persistence ;
    }

    return n ;
}


/////////////////// Fonctions de définition la scène ///////////////////

// Fonction définissant la caméra (code fourni en cours)
Ray camRay (in vec2 p) // p : pixel courant
{
    // Coordonnées normalisées de la souris
    vec2 m = (2. * iMouse.xy - R) / R ;

    // Position de la caméra (exprimée en coordonnées sphériques)
    float theta = m.x * PI ;
    float phi = PHI_DEG * PI / 180. ;
    vec3 ro = CP * vec3(sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta));

    // Axes liés à la caméra
    vec3 ta = vec3(0.0) ;               // point cible
    vec3 cw = normalize(ta-ro) ;        // vecteur "view"
    vec3 cp = vec3(0.0,1.0,0.0) ;       // vecteur "up"
    vec3 cu = normalize(cross(cw,cp)) ; // vecteur "right"
    vec3 cv = normalize(cross(cu,cw)) ; // vecteur "up" normalisé

    // Champ de vue
    float fovDeg = 70. ;
    float fovRad = (fovDeg/180.)*PI;
    float zf = 1./tan(fovRad/2.);
    
    // Vecteur direction normalisé
    vec3 rd = normalize(p.x*cu + p.y*cv*(R.y/R.x) + zf*cw) ;

    return Ray(ro,rd) ;
}

// Fonction renvoyant l'astre correspondant à une planète ou un satellite de la scène
Astre calculAstre (DonneesGeoAstre d, vec3 positionAstreRef, vec3 couleur, vec3 p)
{
    vec3 position = positionAstreRef + d.rayTraj * vec3 (cos(2.*PI*d.freqRev*iTime), 0., sin(2.*PI*d.freqRev*iTime)) ;
    float dist = sdSphere (p + position, d.rayAstre) ;
    return Astre (position, dist, couleur, false) ;
}

// Fonction de la scène 
Astre scene (in vec3 p)
{
    // Tableau des surfaces
    Astre astres[NB_ASTRES] ;

    // Soleil ;
    Astre Soleil = Astre (vec3(0.0), sdSphere(p, RAYON_SOLEIL), vec3(0.), true) ;
    astres[0] = Soleil ;
    
    // Mercure
    float rayonMercure = RT * 0.3825 ;
    float rayonTrajMercure = UA * 0.3871 ;
    float freqMercure = FT * 4.151 ;
    DonneesGeoAstre donneesGeoMercure = DonneesGeoAstre (rayonMercure, rayonTrajMercure, freqMercure) ;
    vec3 couleurMercure = vec3 (114./255., 59./255., 5./255.) ;
    astres[1] = calculAstre (donneesGeoMercure, Soleil.position, couleurMercure, p) ;
    
    // Venus
    float rayonVenus = RT * 0.9488 ;
    float rayonTrajVenus = UA * 0.7230 ;
    float freqVenus = FT * 1.626 ;
    DonneesGeoAstre donneesGeoVenus = DonneesGeoAstre (rayonVenus, rayonTrajVenus, freqVenus) ;
    vec3 couleurVenus = vec3 (226./255., 189./255., 124./255.) ;
    astres[2] = calculAstre (donneesGeoVenus, Soleil.position, couleurVenus, p) ;
      
    // Terre
    float rayonTerre = RT ;
    float rayonTrajTerre = UA ;
    float freqTerre = FT ;
    DonneesGeoAstre donneesGeoTerre = DonneesGeoAstre (rayonTerre, rayonTrajTerre, freqTerre) ;
    vec3 couleurTerre = vec3 (52./255., 61./255., 90./255.) ;
    Astre Terre = calculAstre (donneesGeoTerre, Soleil.position, couleurTerre, p) ;
    astres[3] = Terre ;
    
    // Lune
    float rayonLune = RT * 0.2723 ;
    float rayonTrajLune = RT * 1.5 ;
    float freqLune = FT * 13.37 ;
    DonneesGeoAstre donneesGeoLune = DonneesGeoAstre (rayonLune, rayonTrajLune, freqLune) ;
    vec3 couleurLune = vec3(0.6) ;
    astres[4] = calculAstre (donneesGeoLune, Terre.position, couleurLune, p) ;

    // Mars
    float rayonMars = RT * 0.5325 ;
    float rayonTrajMars = UA * 1.524 ;
    float freqMars = FT * 0.532 ;
    DonneesGeoAstre donneesGeoMars = DonneesGeoAstre (rayonMars, rayonTrajMars, freqMars) ;
    vec3 couleurMars = vec3 (247./255., 129./255., 91./255.) ;
    Astre Mars = calculAstre (donneesGeoMars, Soleil.position, couleurMars, p) ;

    // Cratères de Mars (obtenus par soustraction de sphères qui affleurent la planète)
    int I = 3 ;
    int J = 4 ;
    // I * J cratères tracés
    for (int i=0 ; i<I ; i++)
    {
        for (int j=0 ; j<J ; j++)
        {
            // Calcul d'angles theta et phi aléatoires répartis respectivement entre -PI et PI, puis entre 0 et PI, aux aléa près
            float thetaS = 2.*PI/float(J) * float(j) - PI + 0.25*hash11(float(i*j)) ;
            float phiS = PI * float(i)/float(I) + 0.2*hash11(float(i*j)) ;

            // Ajout à theta d'un angle variant avec le temps pour simuler une rotation selon Y
            float freqRot = 0.5 ;
            float alpha = 2. * PI * freqRot * iTime ;
            thetaS = thetaS + alpha ;

            // Rayon de la sphère et distance du centre de la sphère à Mars (rS)
            float rayonSphere = RT*0.3 + RT*0.1*hash11(float(i*j)) ;
            float rS = rayonMars + rayonSphere - RT * 0.1*hash11(float(i*j*2)) ;
            
            // Position du centre de la sphère en coordonnées sphériques
            vec3 posCentreSphere = rS * vec3 (sin(phiS)*sin(thetaS), cos(phiS), sin(phiS)*cos(thetaS)) ;
            
            // Distance signée à la sphère à soustraire
            float distSphere2 = sdSphere(p+posCentreSphere+Mars.position, rayonSphere) ;

            // Soustraction de la sphère à Mars
            Mars = opDifference(distSphere2, Mars) ;
         }
    }
    astres[5] = Mars ;

    // Jupiter
    float rayonJupiter = RT * 11.21 ;
    float rayonTrajJupiter = UA * 5.203 ;
    float freqJupiter = FT * 0.08430 ;
    DonneesGeoAstre donneesGeoJupiter = DonneesGeoAstre (rayonJupiter, rayonTrajJupiter, freqJupiter) ;
    vec3 couleurJupiter = vec3 (182./255., 96./255., 48./255.) ;
    astres[6] = calculAstre (donneesGeoJupiter, Soleil.position, couleurJupiter, p) ;


    //// Saturne

    // Planète
    float rayonSaturne = RT * 9.449 ;
    float rayonTrajSaturne = UA * 9.582 ;
    float freqSaturne = FT * 0.03395 ;
    DonneesGeoAstre donneesGeoSaturne = DonneesGeoAstre (rayonSaturne, rayonTrajSaturne, freqSaturne) ;
    vec3 couleurSaturne = vec3 (255./255., 212./255., 140./255.) ;
    Astre PlaneteSaturne = calculAstre (donneesGeoSaturne, Soleil.position, couleurSaturne, p) ;
    
    // Anneaux
    vec3 axeSaturne = vec3 (0.05, 0.1, 0.05) ; // axe choisi empiriquement
    vec3 a = -PlaneteSaturne.position ;
    vec3 b = a + axeSaturne ;
      // Anneaux 1
    float rayonMinAnneaux1 = RT * 14.42 ;
    float rayonMaxAnneaux1 = RT * 18.44 ;
    Astre disque1 = Astre (a, sdCylinder(p, a, b, rayonMaxAnneaux1), 0.4*couleurSaturne, false) ;
    Astre sphereSoustraite1 = Astre (a, sdSphere(p-a, rayonMinAnneaux1), vec3(0.), false) ;
    Astre Anneaux1 = opDifference (sphereSoustraite1, disque1) ;
      // Anneaux 2
    float rayonMinAnneaux2 = RT * 19.15 ;
    float rayonMaxAnneaux2 = RT * 21.45 ;
    Astre disque2 = Astre (a, sdCylinder(p, a, b, rayonMaxAnneaux2), 0.25*couleurSaturne, false) ;
    Astre sphereSoustraite2 = Astre (a, sdSphere(p-a, rayonMinAnneaux2), vec3(0.), false) ;
    Astre Anneaux2 = opDifference (sphereSoustraite2, disque2) ;
       // Anneaux complets
    Astre AnneauxSaturne = opUnion(Anneaux1, Anneaux2) ;
    
    // Ensemble planète - anneaux
    astres[7] = opUnion (PlaneteSaturne, AnneauxSaturne) ;


    // Uranus
    float rayonUranus = RT * 4.007 ;
    float rayonTrajUranus = UA * 19.20 ;
    float freqUranus = FT * 0.01190 ;
    DonneesGeoAstre donneesGeoUranus = DonneesGeoAstre (rayonUranus, rayonTrajUranus, freqUranus) ;
    vec3 couleurUranus = vec3 (112./255., 205./255., 254./255.) ;
    astres[8] = calculAstre (donneesGeoUranus, Soleil.position, couleurUranus, p) ;

    // Neptune
    float rayonNeptune = RT * 3.883 ;
    float rayonTrajNeptune = UA * 30.07 ;
    float freqNeptune = FT * 0.006061 ;
    vec3 couleurNeptune = vec3 (111./255., 149./255., 253./255.) ;
    DonneesGeoAstre donneesGeoNeptune = DonneesGeoAstre (rayonNeptune, rayonTrajNeptune, freqNeptune) ;
    astres[9] = calculAstre (donneesGeoNeptune, Soleil.position, couleurNeptune, p) ;

    // Astre "minimal"
    return calculAstreMin(astres) ;
}

// Algorithme de Ray Marching (code fourni en cours)
Astre march (in Ray r, in float distMax)
{
    float d = DIST_MIN ;

    for (int i=0 ; i<RAY_MARCH_STEPS, d<distMax ; ++i)
    {
        Astre a = scene (r.o + d*r.d) ;
        if (a.dist < RAY_MARCH_PRECI)
        {
            return Astre (a.position, d+a.dist, a.couleur, a.type) ;
        }
        d = d + a.dist ;
    }

    return Astre(vec3(0.), distMax, vec3(0), false);
}


/////////// Fonctions liées à la couleur et la lumière ///////////

// Fonction donnant la normale à une surface (code fourni en cours)
vec3 normalAt (in Astre a, in Ray r)
{
    const float e = 0.01 ;
    vec3 p = r.o + a.dist * r.d ;
    float nx = scene(vec3(p.x+e,p.y,p.z)).dist - scene(vec3(p.x-e,p.y,p.z)).dist ;
    float ny = scene(vec3(p.x,p.y+e,p.z)).dist - scene(vec3(p.x,p.y-e,p.z)).dist ;
    float nz = scene(vec3(p.x,p.y,p.z+e)).dist - scene(vec3(p.x,p.y,p.z-e)).dist ;

    return normalize(vec3(nx,ny,nz)) ;
}

// Fonction renvoyant un nombre pseudo-aléatoire pour la variation de l'intensité des étoiles de fond
float scintillement(in float t)
{
    // "Calibration" du scintillement -> a : valeur min, b : valeur max, tau : variation temporelle
    float a = 0.25 ;
    float b = 1. ;
    float tau = 0.5 ;

    // Calculs préalables
    float pt = floor(t) ;
    float f = fract(t) ;
    f = f*f*(3.0-2.0*f) ;

    // Valeur aléatoire de retour
    return a + (b-a) * mix(hash11(pt),hash11(pt+tau),f) ;
}

// Fonction de texture du Soleil
vec3 textureSoleil (in vec3 p)
{
    // Couleurs
    vec3 orange = vec3(1., 140./255., 0.) ;
    vec3 jaune = vec3(1., 1., 0.) ;

    // Bruit
    vec3 t = 40.*p + 2.*iTime ;
    float d = fnoise (t, 0.5, 1., 0.8, 3) ;
    
    return mix (orange, jaune, d) ;
}

// Fonction de shading
vec3 shade (in Astre a, in Ray r)
{
    // Normale au point de la surface atteinte
    vec3 n = normalAt (a,r) ;
    
    // Direction de l'éclairage (le soleil)
    vec3 l = -normalize(r.o + a.dist*r.d);

    // Terme de couleur ambiante
    vec3 Ka = a.couleur * COEF_AMB ;

    // Terme diffus
    vec3 Kd = a.couleur * (1.-COEF_AMB) ;
    float diff = max(dot(n,l),0.) ;

    // Illumination complète (modèle de Phong modifié)
    return Ka + Kd * diff ;
}


/////////////////// Fonction principale ///////////////////

void mainImage (out vec4 fragColor, in vec2 fragCoord)
{
    // Couleurs
    vec3 noir  = vec3(0.) ;
    vec3 blanc = vec3(1.) ;
    
    // Normalisation des coordonnées
    vec2 uv = 2.*(fragCoord / R.xy) - 1. ;

    // Ray marching
    Ray r = camRay (uv) ;
    Astre a = march (r, DIST_MAX) ;
      
    // Couleur du pixel
    vec3 c = noir ;
    if (a.dist < DIST_MAX) // un objet a été atteint
    {
        if (a.type) {  // il s'agit du soleil
            c = textureSoleil (r.o + a.dist*r.d) ;
        }
        else {         // il s'agit de tout autre type d'astre
            // Couleur obtenue par shading
            c = shade(a, r) ;

            // Ombre éventuelle
            vec3 rPointAtteint = r.o + a.dist * r.d ; 
            vec3 n = normalAt(a,r) ;
            vec3 l = -normalize(rPointAtteint) ;
            if (dot(n,l) > 0.0)  // on ajoute de l'ombre seulement pour les points éclairés par le soleil
            {
                float distMax = length(rPointAtteint)-(RAYON_SOLEIL+0.001) ; // on recherche des objets entre l'astre atteint et le soleil, pas plus loin 
                Ray r2 = Ray(rPointAtteint, l) ;
                Astre a2 = march(r2, distMax) ;
                if (a2.dist < distMax) {
                    c = a.couleur * COEF_AMB ;
                }
            }
        }
    }
    else  // pas d'objet atteint, on ajoute (éventuellement) une étoile
    {
        float seuil = 0.0001 ;
        float decalageAlea = 1000. * hash11(uv.x+uv.y) ;
        c = hash12(uv) < seuil ? scintillement(iTime+decalageAlea)*blanc : noir ;
    }
    
    fragColor = vec4(c,1.0) ;
}


