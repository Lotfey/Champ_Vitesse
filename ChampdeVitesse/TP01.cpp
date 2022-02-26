// TP01.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "../CImg.h"

using namespace std;
using namespace cimg_library;


CImg<double> HornSchunck(CImg<unsigned char> seq)
{
 CImg<double> field(seq.dimx(),seq.dimy(),1,2);
 CImg<double> field_conv;

 // Creation du noyau de filtrage pour le calcul du Laplacien
 CImg<double> Mt(3,3); 
 Mt(0,0) = .25;	Mt(0,1) = .5;	Mt(0,2)= .25;
 Mt(1,0) = .5;	Mt(1,1) = 0.;	Mt(1,2)= .5;
 Mt(2,0) = .25;	Mt(2,1) = .5;	Mt(2,2)= .25;

 // Initialisation du champ
 field.fill(0);

 // Calcul des gradients de l'image
 CImgList<double> grad = seq.get_gradient("xyz", +1);
 
 // Schéma iteratif
 int niter = 10;
 double alpha = 100, alpha2 = alpha * alpha;
 for (int n = 0; n < niter; n++)
 {
	CImg<double> moy = field.get_convolve(Mt);

	for (int x = 0; x < field.dimx() - 1; x++)
		for (int y = 0; y < field.dimy() - 1; y++)
		{
			double Ix = 0.25 * (grad[0](x,y,0) + grad[0](x,y+1,0) + grad[0](x,y,1) + grad[0](x,y+1,1));
			double Iy = 0.25 * (grad[1](x,y,0) + grad[1](x+1,y,0) + grad[1](x,y,1) + grad[1](x+1,y,1));
			double It = 0.25 * (grad[2](x,y,0) + grad[2](x+1,y,0) + grad[2](x,y+1,0) + grad[2](x+1,y+1,0));
			double un = moy(x, y, 0), vn = moy(x, y, 1);
			double num = Ix * un + Iy * vn + It;
			double den = alpha2 + Ix * Ix + Iy * Iy;
			field(x, y, 0, 0) = un - Ix * num / den;
			field(x, y, 0, 1) = vn - Iy * num / den;
		}
 }

 return field;
}

CImg<double> LucasKanade(CImg<unsigned char> seq)
{
 int n = 9;				   // Taille du voisinage
 CImg<double> A(2,2); 
 CImg<double> b(2);   
 CImg<double> weight(n,n); 

 // Allocation du champ de sortie
 CImg<double> field(seq.dimx(),seq.dimy(),1,2);
 field.fill(0);

 // Calcul des gradients
 CImgList<double> grad = seq.get_gradient("xyz",4);

 // Calcul de la matrice de ponderation
 float sigma = 10;
 float color = 1;
 weight.draw_gaussian((float)n/2,(float)n/2,sigma,&color);

 // Calcul du champ en chaque points
 for (int x = n / 2; x < field.dimx() - n / 2; x++)
	 for (int y = n / 2; y < field.dimy() - n / 2; y++)
	 {
		 A.fill(0.0);
		 b.fill(0.0);
		 for (int i = x - n / 2; i <= x + n / 2; i++)
			 for (int j = y - n / 2; j <= y + n / 2; j++)
			 {
				 double Ix = grad[0](i, j, 0);
				 double Iy = grad[1](i, j, 0);
				 double It = grad[2](i, j, 0);
				 A(0, 0) += Ix * Ix;
				 A(0, 1) += Ix * Iy;
				 A(1, 1) += Iy * Iy;
				 b(0) -= Ix * It;
				 b(1) -= Iy * It;
			 }
		 A(1, 0) = A(0, 1);

		 // Résolution du système 2x2
		 double det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
		 if (fabs(det) > 1.e-6)
		 {
			 double u = (b(0) * A(1, 1) - b(1) * A(0, 1)) / det;
			 double v = (b(0) - A(0, 0) * u) / A(0, 1);
			 field(x, y, 0, 0) = u;
			 field(x, y, 0, 1) = v;
		 }
	 }

 
 return field;
}

int _tmain(int argc, _TCHAR* argv[])
{
 // Ouverture des images d'entrée
 CImg<unsigned char> seq = CImg<>("../Img/4_cav_1.bmp").channel(0);
 seq.append(CImg<unsigned char>("../Img/4_cav_2.bmp").channel(0),'z');

 CImg<> disp;
 char str[256];
 
 int choice;
 cout << "Choix de l'algorithme (0:Horn et Schunk, 1:Lucas et Kanade) : ";
 cin >> choice;

 // Horn et Schunk
 if (choice == 0)
 {
	strcpy(str,"Horn et Schunk");
	disp = HornSchunck(seq);
 }
 
 // Lucas et Kanade
 if (choice == 1)
 {
	strcpy(str,"Lucas et Kanade");
	disp = LucasKanade(seq);
 }

 // Affichage du champ résultat
 float color=500; unsigned int  sampling = 8; float factor = 5; int  quiver_type = 0; float  opacity = 0.5;

 CImg<unsigned char> out1 = seq.get_slice(0).draw_quiver(disp,&color,opacity,sampling,factor,quiver_type);
 out1.save("lucas_4cav.bmp");
 CImgDisplay disp_res1(out1,str);
 CImg<unsigned char> out2 = seq.get_slice(0);
 CImgDisplay disp_res2(out2,"Image t");
 CImg<unsigned char> out3 = seq.get_slice(1);
 CImgDisplay disp_res3(out3,"Image t+1");

 while (!disp_res1.is_closed ) disp_res1.wait();

 return 0;
}
