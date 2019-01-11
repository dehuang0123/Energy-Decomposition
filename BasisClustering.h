#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include "PairingPackage.h"


void Clustering(std::vector<Vertex> *Vertex_list, int N, int M)
{
	std::vector<Patch> Patch_list(N);
	std::vector<Pedge> Pedge_list(M);
	std::vector<Patch *> Patch_key(N);
	std::vector<PVedge> Patch_vertex_list(N);
	int i, j, l, k = 0;
	int tempDegree;
	int eindex;
	double selfcoe;
	Vedge *Ve_temp;
	VEedge *VE_temp;
	Pedge *Pe_temp;
	PVedge *PVe_temp;


	for (i = 0; i<N; ++i)
	{
		Patch_list[i].Index = i;
	}

	// Initialize Patch_list
	for (i = 0; i<N; ++i)
	{
		Patch_key[i] = &Patch_list[i];

		
		Patch_list[i].Mark = 0;
		Patch_list[i].Cluster = &Patch_list[i];

		// Initialize members
		Patch_list[i].Psize = 1;
		Patch_vertex_list[i].Outvertex = &(*Vertex_list)[i];
		Patch_vertex_list[i].Next = NULL;
		Patch_list[i].Member = &Patch_vertex_list[i];
		Patch_list[i].Lastmember = &Patch_vertex_list[i];

		// Initialize neighbours
		Patch_list[i].Neighbour = &Pedge_list[k];
		Patch_list[i].Delta = 0;
		VE_temp = (*Vertex_list)[i].Interact; 
		tempDegree = 0;

		for (j = 0; j<(*Vertex_list)[i].Edegree; ++j)
		{
			eindex = VE_temp->Eindex;

			while ( VE_temp != NULL && VE_temp-> Eindex == eindex )
			{

				if ( VE_temp->Outvertex->Index == i )
				{
					Patch_list[i].Delta += std::abs(VE_temp-> Coefficient);
				}
				else 
				{

					for ( l = 0; l<tempDegree ; ++l )
					{

						if ( VE_temp->Outvertex->Index == Pedge_list[k+l].Outpatch->Index )
						{
							Pedge_list[k+l].Bond += std::abs(VE_temp-> Coefficient);
							Patch_list[i].Delta += std::abs(VE_temp-> Coefficient);
							break;
						}
					}

					if ( l == tempDegree )
					{

						tempDegree += 1;
						Pedge_list[k+l].Outpatch = &Patch_list[VE_temp->Outvertex->Index];
						Pedge_list[k+l].Bond = std::abs(VE_temp-> Coefficient);
						Patch_list[i].Delta += std::abs(VE_temp-> Coefficient);
						if ( l > 0 ) 
						{
							Pedge_list[k+l-1].Next = &Pedge_list[k+l];
						}
					}
				}

				VE_temp = VE_temp->Next;

			}
		}

		k += tempDegree; 

		Patch_list[i].Degree = tempDegree;
		Patch_list[i].Lastneighbour = &Pedge_list[k - 1];
	}

	std::cout << N << " " << M << std::endl;

	int Patch_N = 0, s = 0, Pedge_M = 0;

	Patch_N = Parallel_Pairing(&Patch_key, N);

	std::cout << std::endl;
	std::cout << Patch_N << std::endl;

	std::ofstream ClusterOutput;
	std::ofstream NeighbourOutput;

	ClusterOutput.open("cluster.txt");
	NeighbourOutput.open("neighbour.txt");

	ClusterOutput << Patch_N << "\n";
	NeighbourOutput << Patch_N << "\n";

	double max_delta = 0;

	for (i = 0; i<Patch_N; ++i)
	{
		Pedge_M += Patch_key[i]->Degree;
		s += Patch_key[i]->Psize;
		ClusterOutput << Patch_key[i]->Index << " ";
		ClusterOutput << Patch_key[i]->Delta << " ";

		NeighbourOutput << Patch_key[i]->Index << " ";

		if (max_delta<  Patch_key[i]->Delta) max_delta = Patch_key[i]->Delta;

		ClusterOutput << Patch_key[i]->Lambda << " ";
		ClusterOutput << "\n";
		ClusterOutput << Patch_key[i]->Psize << " ";
		PVe_temp = Patch_key[i]->Member;
		for (j = 0; j<Patch_key[i]->Psize; ++j)
		{
			ClusterOutput << PVe_temp->Outvertex->Index << " ";
			ClusterOutput << PVe_temp->Outvertex->Eigfun << " ";
			PVe_temp = PVe_temp->Next;
		}
		ClusterOutput << "\n";

		NeighbourOutput<< Patch_key[i]->Degree <<"\n";
		Pe_temp=Patch_key[i]->Neighbour;
		for (j=0;j<Patch_key[i]->Degree;++j)
		{
			NeighbourOutput<<Pe_temp->Outpatch->Index<<" ";
			Pe_temp=Pe_temp->Next;
		} 		
		NeighbourOutput<<"\n";
	}

	ClusterOutput.close();
	NeighbourOutput.close();

	std::cout << max_delta << std::endl;

	// for (i=0;i<Patch_N;++i)
	// {	
	// 	s=s+Patch_key[i]->Psize;
	// 	std::cout<< Patch_key[i]->Psize <<" ";
	// 	Ve_temp=Patch_key[i]->Member;
	// for (j=0;j<Patch_key[i]->Psize;++j)
	// {
	// 	std::cout<<Ve_temp->Outvertex->Index<<" ";
	// 	Ve_temp=Ve_temp->Next;
	// }
	// std::cout<<std::endl;
	// }
	std::cout << s << std::endl;

	// for (i=0;i<N;++i)
	// {
	// 	std::cout<< Patch_key[i]->Index<<" "<< Patch_list[i].Delta<< " ";
	// 	Pe_temp=Patch_list[i].Neighbour;
	// 	for (j=0;j<Vertex_list[i].Degree;++j)
	// 	{
	// 		std::cout<< Pe_temp->Bond<< " ";
	// 		Pe_temp=Pe_temp->Next;
	// 	}
	// 	std::cout<< std::endl;
	// }

}