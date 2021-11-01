#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

/*
 * ARGS TO PASS:
 * ~~~~~~~~~~~~~
 *
 * argv[0] = program
 * argv[1] = input dump file name
 * argv[2] = number of water molecules (these molecules will be added randomly)
 *
 */

int getNatoms (const char *inputFileName)
{
	FILE *read = fopen (inputFileName, "r");
	char lineString[2000];
	int lineNumber = 0, natoms;

	while (fgets (lineString, 2000, read) != NULL)
	{
		lineNumber++;

		if (lineNumber == 4)
		{
			sscanf (lineString, "%d", &natoms);
			return natoms;
		}
	}
	return 0;
}

typedef struct datafile_atoms
{
	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct dump
{
	int id, type, ix, iy, iz;
	float x, y, z, xs, ys, zs;
} DUMP;

typedef struct bounds
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
} BOUNDS;

float findMaxChainDimension (const char *inputFile, int nAtoms)
{
	FILE *input;
	input = fopen (inputFile, "r");
	int nTimeframes = 0, arrayIndex = 0;
	float maxDimension = 0, maxDimension_avg = 0, dimension;
	char lineString[2000];

	DUMP *traj, *traj_temp, com, dumpLow, dumpHigh;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	traj_temp = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	com.x = 0; com.y = 0; com.z = 0;

	// Get simulation box dimensions from dump file
	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 1000, input);
	}
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.x, &dumpHigh.x);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.y, &dumpHigh.y);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.z, &dumpHigh.z);
	rewind (input);

	printf("Computing max chain dimension...\n");

	while ((fgets (lineString, 2000, input) != NULL))
	{
		maxDimension = 0;
		if (strstr (lineString, "ITEM: ATOMS id type x y z"))
		{
			for (int i = 0; i < nAtoms; ++i)
			{
				fgets (lineString, 2000, input);
				sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", &traj[i].id, &traj[i].type, &traj[i].x, &traj[i].y, &traj[i].z, &traj[i].xs, &traj[i].ys, &traj[i].zs, &traj[i].ix, &traj[i].iy, &traj[i].iz);

				traj[i].x = traj[i].x + traj[i].ix * (dumpHigh.x - dumpLow.x);
				traj[i].y = traj[i].y + traj[i].iy * (dumpHigh.y - dumpLow.y);
				traj[i].z = traj[i].z + traj[i].iz * (dumpHigh.z - dumpLow.z);

				com.x += traj[i].x;
				com.y += traj[i].y;
				com.z += traj[i].z;
			}

			com.x /= nAtoms; com.y /= nAtoms; com.z /= nAtoms;

			for (int i = 0; i < nAtoms; ++i)
			{
				dimension = sqrt (pow ((traj[i].x - com.x), 2) + pow ((traj[i].y - com.y), 2) + pow ((traj[i].z - com.z), 2));
				if (dimension > maxDimension)
					maxDimension = dimension;
			}
			maxDimension_avg += maxDimension;

			nTimeframes++;
			if (!(nTimeframes%10))
			{
				fprintf(stdout, "\rScanning timeframe: %d", nTimeframes);
				fflush (stdout);
			}
		}
	}

	printf("\n\n");

	maxDimension_avg /= nTimeframes;

	// Previous calculation only gives the maxDimension from the center of mass
	// So, the calculated value is multiplied by 2 before returning to the main function
	maxDimension_avg *= 2;

	printf("\nmaxDimension: %f\n", maxDimension_avg);

	return maxDimension_avg;
}

DUMP findMaxChainDimension_XYZ (const char *inputFile, int nAtoms)
{
	FILE *input;
	input = fopen (inputFile, "r");
	int nTimeframes = 0, arrayIndex = 0;
	float maxDimension_x = 0, maxDimension_y = 0, maxDimension_z = 0, maxDimension_avg_x = 0, maxDimension_avg_y = 0, maxDimension_avg_z = 0, dimension_x, dimension_y, dimension_z;
	char lineString[2000];

	DUMP *traj, *traj_temp, com, dumpLow, dumpHigh, chainDimension;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	traj_temp = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	com.x = 0; com.y = 0; com.z = 0;

	// Get simulation box dimensions from dump file
	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 1000, input);
	}
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.x, &dumpHigh.x);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.y, &dumpHigh.y);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpLow.z, &dumpHigh.z);
	rewind (input);

	printf("Computing max chain dimension...\n");

	while ((fgets (lineString, 2000, input) != NULL))
	{
		maxDimension_x = 0;
		maxDimension_y = 0;
		maxDimension_z = 0;

		if (strstr (lineString, "ITEM: ATOMS id type x y z"))
		{
			for (int i = 0; i < nAtoms; ++i)
			{
				fgets (lineString, 2000, input);
				sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", &traj[i].id, &traj[i].type, &traj[i].x, &traj[i].y, &traj[i].z, &traj[i].xs, &traj[i].ys, &traj[i].zs, &traj[i].ix, &traj[i].iy, &traj[i].iz);

				traj[i].x = traj[i].x + traj[i].ix * (dumpHigh.x - dumpLow.x);
				traj[i].y = traj[i].y + traj[i].iy * (dumpHigh.y - dumpLow.y);
				traj[i].z = traj[i].z + traj[i].iz * (dumpHigh.z - dumpLow.z);

				com.x += traj[i].x;
				com.y += traj[i].y;
				com.z += traj[i].z;
			}

			com.x /= nAtoms; com.y /= nAtoms; com.z /= nAtoms;

			for (int i = 0; i < nAtoms; ++i)
			{
				dimension_x = sqrt (pow ((traj[i].x - com.x), 2));
				dimension_y = sqrt (pow ((traj[i].y - com.y), 2));
				dimension_z = sqrt (pow ((traj[i].z - com.z), 2));
				if (dimension_x > maxDimension_x)
					maxDimension_x = dimension_x;
				if (dimension_y > maxDimension_y)
					maxDimension_y = dimension_y;
				if (dimension_z > maxDimension_z)
					maxDimension_z = dimension_z;
			}
			maxDimension_avg_x += maxDimension_x;
			maxDimension_avg_y += maxDimension_y;
			maxDimension_avg_z += maxDimension_z;

			nTimeframes++;
			if (!(nTimeframes%10))
			{
				fprintf(stdout, "\rScanning timeframe: %d", nTimeframes);
				fflush (stdout);
			}
		}
	}

	printf("\n\n");

	maxDimension_avg_x /= nTimeframes;
	maxDimension_avg_y /= nTimeframes;
	maxDimension_avg_z /= nTimeframes;

	// Previous calculation only gives the maxDimension from the center of mass
	// So, the calculated value is multiplied by 2 before returning to the main function

	chainDimension.x = maxDimension_avg_x * 2;
	chainDimension.y = maxDimension_avg_y * 2;
	chainDimension.z = maxDimension_avg_z * 2;

	printf(" maxDimension_avg_x: %f\n maxDimension_avg_y: %f\n maxDimension_avg_z: %f\n", chainDimension.x, chainDimension.y, chainDimension.z);

	return chainDimension;
}

DUMP findCOM (DUMP *traj, int lineCount)
{
	DUMP com;
	com.x = 0; com.y = 0; com.z = 0;

	for (int i = 0; i < lineCount; ++i)
	{
		com.x += traj[i].x;
		com.y += traj[i].y;
		com.z += traj[i].z;
	}

	com.x /= lineCount;
	com.y /= lineCount;
	com.z /= lineCount;

	return com;
}

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	FILE *input;
	input = fopen (dataFileName, "r");

	int isAtomLine = 0, nAtoms = -1, nAtomLine = 0;
	int isBondLine = 0, nBonds = -1, nBondLine = 0;
	int isAngleLine = 0, nAngles = -1, nAngleLine = 0;
	int isDihedralLine = 0, nDihedrals = -1, nDihedralLine = 0;
	int isImproperLine = 0, nImpropers = -1, nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	// DATA_ATOMS *atoms;
	// DATA_BONDS *bonds;
	// DATA_ANGLES *angles;
	// DATA_DIHEDRALS *dihedrals;
	// DATA_IMPROPERS *impropers;
	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);

		if (strstr (lineString, "bond types"))
			sscanf (lineString, "%d \n", &datafile.nBondTypes);

		if (strstr (lineString, "angle types"))
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);

		if (strstr (lineString, "dihedral types"))
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);

		if (strstr (lineString, "improper types"))
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	return datafile;
}

DUMP *readLastFrame (const char *pipeString, int nAtoms, BOUNDS dumpDimension)
{
	FILE *input;
	input = popen (pipeString, "r");

	DUMP *traj, *traj_temp, com, max;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	traj_temp = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	int lineCount = 0;

	char lineString[1000];

	while (fgets (lineString, 1000, input) != NULL)
	{
		sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", &traj_temp[lineCount].id, &traj_temp[lineCount].type, &traj_temp[lineCount].x, &traj_temp[lineCount].y, &traj_temp[lineCount].z, &traj_temp[lineCount].xs, &traj_temp[lineCount].ys, &traj_temp[lineCount].zs, &traj_temp[lineCount].ix, &traj_temp[lineCount].iy, &traj_temp[lineCount].iz);
		lineCount++;
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nAtoms; ++j)
		{
			if (traj_temp[j].id == i + 1)
			{
				traj[i].id = traj_temp[j].id;
				traj[i].type = traj_temp[j].type;
				traj[i].x = traj_temp[j].x;
				traj[i].y = traj_temp[j].y;
				traj[i].z = traj_temp[j].z;
				traj[i].xs = traj_temp[j].xs;
				traj[i].ys = traj_temp[j].ys;
				traj[i].zs = traj_temp[j].zs;
				traj[i].ix = traj_temp[j].ix;
				traj[i].iy = traj_temp[j].iy;
				traj[i].iz = traj_temp[j].iz;
			}
		}
	}

	// Finding the highest image flag
	max.ix = 0; max.iy = 0; max.iz = 0;
	for (int i = 0; i < nAtoms; ++i)
	{
		if (abs (traj[i].ix) > abs (max.ix))
			max.ix = traj[i].ix;
		if (abs (traj[i].iy) > abs (max.iy))
			max.iy = traj[i].iy;
		if (abs (traj[i].iz) > abs (max.iz))
			max.iz = traj[i].iz;
	}

	printf("\nSimulation box image from trajectory file.\n\n   max.ix: %d\n   max.iy: %d\n   max.iz: %d\n\n", max.ix, max.iy, max.iz);

	// Rewriting image information
	for (int i = 0; i < nAtoms; ++i)
	{
		traj[i].ix -= max.ix;
		traj[i].iy -= max.iy;
		traj[i].iz -= max.iz;
	}

	// Unwrap coordinates
	for (int i = 0; i < nAtoms; ++i)
	{
		if (traj[i].ix > 0)
			traj[i].x = (traj[i].x - dumpDimension.xlo) + dumpDimension.xhi;
		else if (traj[i].ix < 0)
			traj[i].x = dumpDimension.xlo - (dumpDimension.xhi - traj[i].x);
		if (traj[i].iy > 0)
			traj[i].y = (traj[i].y - dumpDimension.ylo) + dumpDimension.yhi;
		else if (traj[i].iy < 0)
			traj[i].y = dumpDimension.ylo - (dumpDimension.yhi - traj[i].y);
		if (traj[i].iz > 0)
			traj[i].z = (traj[i].z - dumpDimension.zlo) + dumpDimension.zhi;
		else if (traj[i].iz < 0)
			traj[i].z = dumpDimension.zlo - (dumpDimension.zhi - traj[i].z);
	}

	// Finding the center of mass
	com.x = 0; com.y = 0; com.z = 0;
	for (int i = 0; i < nAtoms; ++i)
	{
		com.x += traj[i].x;
		com.y += traj[i].y;
		com.z += traj[i].z;
	}

	com.x /= nAtoms; com.y /= nAtoms; com.z /= nAtoms;
	printf("Centers of traj from input dump file (last timeframe)\n com.x: %f\n com.y: %f\n com.z: %f\n", com.x, com.y, com.z);

	// Recentering the chain
	for (int i = 0; i < nAtoms; ++i)
	{
		traj[i].x -= com.x;
		traj[i].y -= com.y;
		traj[i].z -= com.z;
	}

	pclose (input);
	return traj;
}

BOUNDS getDumpDimension (FILE *input)
{
	BOUNDS dumpDimension;

	char lineString[1000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 1000, input);
	}

	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpDimension.xlo, &dumpDimension.xhi);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpDimension.ylo, &dumpDimension.yhi);
	fgets (lineString, 1000, input);
	sscanf (lineString, "%f %f\n", &dumpDimension.zlo, &dumpDimension.zhi);

	printf("Simulation box dimensions from input dump file:\n\n xlo xhi %f %f\n ylo yhi %f %f\n zlo zhi %f %f\n", dumpDimension.xlo, dumpDimension.xhi, dumpDimension.ylo, dumpDimension.yhi, dumpDimension.zlo, dumpDimension.zhi);

	rewind (input);

	return dumpDimension;
}

DUMP *recenterCoords (DUMP *traj_in, DUMP com, BOUNDS dumpDimension, int nAtoms)
{
	DUMP *traj;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	// Recenter chain
	printf("Recentering chain from dump file...\n");
	for (int i = 0; i < nAtoms; ++i)
	{
		traj[i].x = traj_in[i].x + (-1 * com.x);
		traj[i].y = traj_in[i].y + (-1 * com.y);
		traj[i].z = traj_in[i].z + (-1 * com.z);
	}

	return traj;
}

BOUNDS computeSimBoxDimension (DUMP chainDimension)
{
	BOUNDS simBoxDimension;
	float scaleX, scaleY, scaleZ;

	printf("%s\n", "Simulation box dimension will be decided based on max chain dimension according to the expression below...\n xlo = maxDimension * -x\n xhi = maxDimension * x\n ylo = maxDimension * -y\n yhi = maxDimension * y\n zlo = maxDimension * -z\n zhi = maxDimension * z\n");
	printf("%s ", "Enter x:"); scanf ("%f", &scaleX); 
	printf("%s ", "Enter y:"); scanf ("%f", &scaleY); 
	printf("%s ", "Enter z:"); scanf ("%f", &scaleZ); 

	simBoxDimension.xlo = (chainDimension.x / 2) * -scaleX;
	simBoxDimension.xhi = (chainDimension.x / 2) * scaleX;
	simBoxDimension.ylo = (chainDimension.y / 2) * -scaleY;
	simBoxDimension.yhi = (chainDimension.y / 2) * scaleY;
	simBoxDimension.zlo = (chainDimension.z / 2) * -scaleZ;
	simBoxDimension.zhi = (chainDimension.z / 2) * scaleZ;

	printf("\nSimulation box dimension (recalculated)\n\n xlo: %f\txhi: %f (length: %f)\n ylo: %f\tyhi: %f (length: %f)\n zlo: %f\tzhi: %f (length: %f)\n", simBoxDimension.xlo, simBoxDimension.xhi, fabs (simBoxDimension.xhi) + fabs (simBoxDimension.xlo), simBoxDimension.ylo, simBoxDimension.yhi, fabs (simBoxDimension.yhi) + fabs (simBoxDimension.ylo), simBoxDimension.zlo, simBoxDimension.zhi, fabs (simBoxDimension.zhi) + fabs (simBoxDimension.zlo));

	return simBoxDimension;
}

int findNButanediol (float butanediolFraction, BOUNDS simBoxDimension)
{
	float xLength = fabs (simBoxDimension.xhi) + fabs (simBoxDimension.xlo), yLength = fabs (simBoxDimension.yhi) + fabs (simBoxDimension.ylo), zLength = fabs (simBoxDimension.zhi) + fabs (simBoxDimension.zlo), avogadroNumber = 6.023, butanediolDensity = 1.02, molarMass = 90.12;
	float nButanediol_max = (butanediolDensity * xLength * yLength * zLength * avogadroNumber * 0.1 / molarMass);
	int nButanediol = (int) (nButanediol_max * butanediolFraction);
	printf("%d molecules of butanediol will be added... (nButanediol_max: %.0f)\n", nButanediol, nButanediol_max);

	return nButanediol;
}

int findNWater (float waterFraction, BOUNDS simBoxDimension)
{
	printf("Volume fraction of water in the simulation box is set as: %.2f\n", waterFraction);
	float xLength = fabs (simBoxDimension.xhi) + fabs (simBoxDimension.xlo), yLength = fabs (simBoxDimension.yhi) + fabs (simBoxDimension.ylo), zLength = fabs (simBoxDimension.zhi) + fabs (simBoxDimension.zlo), avogadroNumber = 6.023, waterDensity = 1.0, molarMass = 18;
	float nWater_max = (waterDensity * xLength * yLength * zLength * avogadroNumber * 0.1 / molarMass);
	int nWater = (int) (nWater_max * waterFraction);
	printf("%d molecules of water will be added... (nWater_max: %.0f)\n", nWater, nWater_max);

	return nWater;
}

DATA_ATOMS *populateButanediol (BOUNDS simBoxDimension, int nButanediol)
{
	DATA_ATOMS *butanediol;
	butanediol = (DATA_ATOMS *) calloc (nButanediol, sizeof (DATA_ATOMS));

	float nBins_x = cbrt (nButanediol), nBins_y = cbrt (nButanediol), nBins_z =cbrt (nButanediol);
	float xDistSeparation = (simBoxDimension.xhi + simBoxDimension.xlo) / nBins_x, yDistSeparation = (simBoxDimension.yhi + simBoxDimension.ylo) / nBins_y, zDistSeparation = (simBoxDimension.zhi + simBoxDimension.zlo) / nBins_z;
	int currentButanediol = 0;

	// Distribute the butanediol molecules evenly
	for (int i = 0; i < nBins_x; ++i)
	{
		for (int j = 0; j < nBins_y; ++j)
		{
			for (int k = 0; k < nBins_z; ++k)
			{
				butanediol[currentButanediol].x = simBoxDimension.xlo + ((i + 1) * xDistSeparation) - (xDistSeparation / 2);
				butanediol[currentButanediol].y = simBoxDimension.ylo + ((j + 1) * yDistSeparation) - (yDistSeparation / 2);
				butanediol[currentButanediol].z = simBoxDimension.zlo + ((k + 1) * zDistSeparation) - (zDistSeparation / 2);
				currentButanediol++;
			}
		}
	}

	return butanediol;
}

DATA_ATOMS *populateWater (BOUNDS simBoxDimension, int nWater)
{
	DATA_ATOMS *water;
	water = (DATA_ATOMS *) calloc (nWater, sizeof (DATA_ATOMS));

	float nBins_x = cbrt (nWater), nBins_y = cbrt (nWater), nBins_z = cbrt (nWater);
	float xDistSeparation = (fabs (simBoxDimension.xhi) + fabs (simBoxDimension.xlo)) / nBins_x, yDistSeparation = ( fabs(simBoxDimension.yhi) + fabs (simBoxDimension.ylo)) / nBins_y, zDistSeparation = (fabs (simBoxDimension.zhi) + fabs (simBoxDimension.zlo)) / nBins_z;
	int currentWater = 0;

	// Distributing water molecules evenly
	for (int i = 0; i < nBins_x; ++i)
	{
		for (int j = 0; j < nBins_y; ++j)
		{
			for (int k = 0; k < nBins_z; ++k)
			{
				water[currentWater].x = simBoxDimension.xlo + ((i + 1) * xDistSeparation) - (xDistSeparation / 2);
				water[currentWater].y = simBoxDimension.ylo + ((j + 1) * yDistSeparation) - (yDistSeparation / 2);
				water[currentWater].z = simBoxDimension.zlo + ((k + 1) * zDistSeparation) - (zDistSeparation / 2);
				currentWater++;
			}
		}
	}

	return water;
}

void createWaterTopology (DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, int nWater, DATA_ATOMS *water, DATAFILE_INFO datafile_raw, DATAFILE_INFO datafile)
{
	int currentIDAtoms = datafile_raw.nAtoms, currentIDBonds = datafile_raw.nBonds, currentIDAngles = datafile_raw.nAngles;
	int O_ID, H1_ID, H2_ID;

	for (int i = 0; i < nWater; ++i)
	{
		// Adding oxygen molecule
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		O_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 2;
		(*atoms)[currentIDAtoms].atomType = 8;
		(*atoms)[currentIDAtoms].charge = -0.83;
		(*atoms)[currentIDAtoms].x = water[i].x;
		(*atoms)[currentIDAtoms].y = water[i].y;
		(*atoms)[currentIDAtoms].z = water[i].z;
		currentIDAtoms++;

		// Adding hydrogens
		// Hydrogen 1
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		H1_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 2;
		(*atoms)[currentIDAtoms].atomType = 9;
		(*atoms)[currentIDAtoms].charge = 0.415;
		(*atoms)[currentIDAtoms].x = water[i].x - 0.32;
		(*atoms)[currentIDAtoms].y = water[i].y + 0.13;
		(*atoms)[currentIDAtoms].z = water[i].z - 0.91;
		currentIDAtoms++;

		// Hydrogen 2
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		H2_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 2;
		(*atoms)[currentIDAtoms].atomType = 9;
		(*atoms)[currentIDAtoms].charge = 0.415;
		(*atoms)[currentIDAtoms].x = water[i].x + 0.97;
		(*atoms)[currentIDAtoms].y = water[i].y;
		(*atoms)[currentIDAtoms].z = water[i].z;
		currentIDAtoms++;

		// Adding bonds (bond1 and bond2 are identical)
		// Bond 1
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 7;
		(*bonds)[currentIDBonds].atom1 = H2_ID;
		(*bonds)[currentIDBonds].atom2 = O_ID;
		currentIDBonds++;

		// Bond 2
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 7;
		(*bonds)[currentIDBonds].atom1 = H1_ID;
		(*bonds)[currentIDBonds].atom2 = O_ID;
		currentIDBonds++;

		// Adding angles
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 10;
		(*angles)[currentIDAngles].atom1 = H2_ID;
		(*angles)[currentIDAngles].atom2 = O_ID;
		(*angles)[currentIDAngles].atom3 = H1_ID;
		currentIDAngles++;
	}
}

void createButanediolTopology (DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, int nButanediol, DATA_ATOMS *butanediol, DATAFILE_INFO datafile_raw, DATAFILE_INFO datafile)
{
	int currentIDAtoms = datafile.nAtoms, currentIDBonds = datafile.nBonds, currentIDAngles = datafile.nAngles, currentIDDihedrals = datafile.nDihedrals;

	int H1_ID, O1_ID, C1_ID, C2_ID, C3_ID, C4_ID, O2_ID, H2_ID;

	for (int i = 0; i < nButanediol; ++i)
	{
		// Defined by TraPPE force field
		// Adding atomic coordinates
		// H
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		H1_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 10;
		(*atoms)[currentIDAtoms].charge = 0.435;
		(*atoms)[currentIDAtoms].x = butanediol[i].x - 2.18;
		(*atoms)[currentIDAtoms].y = butanediol[i].y + 1.29;
		(*atoms)[currentIDAtoms].z = butanediol[i].z - 1.26;
		currentIDAtoms++;

		// O
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		O1_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 11;
		(*atoms)[currentIDAtoms].charge = -0.7;
		(*atoms)[currentIDAtoms].x = butanediol[i].x - 1.93;
		(*atoms)[currentIDAtoms].y = butanediol[i].y + 0.34;
		(*atoms)[currentIDAtoms].z = butanediol[i].z - 1.36;
		currentIDAtoms++;

		// CH2
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		C1_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 12;
		(*atoms)[currentIDAtoms].charge = 0.265;
		(*atoms)[currentIDAtoms].x = butanediol[i].x - 0.53;
		(*atoms)[currentIDAtoms].y = butanediol[i].y + 0.34;
		(*atoms)[currentIDAtoms].z = butanediol[i].z - 1.37;
		currentIDAtoms++;

		// CH2
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		C2_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 12;
		(*atoms)[currentIDAtoms].charge = 0.0;
		(*atoms)[currentIDAtoms].x = butanediol[i].x;
		(*atoms)[currentIDAtoms].y = butanediol[i].y;
		(*atoms)[currentIDAtoms].z = butanediol[i].z;
		currentIDAtoms++;

		// CH2
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		C3_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 12;
		(*atoms)[currentIDAtoms].charge = 0.0;
		(*atoms)[currentIDAtoms].x = butanediol[i].x + 1.52;
		(*atoms)[currentIDAtoms].y = butanediol[i].y;
		(*atoms)[currentIDAtoms].z = butanediol[i].z;
		currentIDAtoms++;

		// CH2
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		C4_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 12;
		(*atoms)[currentIDAtoms].charge = 0.265;
		(*atoms)[currentIDAtoms].x = butanediol[i].x + 2.05;
		(*atoms)[currentIDAtoms].y = butanediol[i].y - 1.25;
		(*atoms)[currentIDAtoms].z = butanediol[i].z - 0.68;
		currentIDAtoms++;

		// O
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		O2_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 11;
		(*atoms)[currentIDAtoms].charge = -0.7;
		(*atoms)[currentIDAtoms].x = butanediol[i].x + 3.44;
		(*atoms)[currentIDAtoms].y = butanediol[i].y - 1.23;
		(*atoms)[currentIDAtoms].z = butanediol[i].z - 0.67;
		currentIDAtoms++;

		// H
		(*atoms)[currentIDAtoms].id = (*atoms)[currentIDAtoms-1].id + 1;
		H2_ID = (*atoms)[currentIDAtoms-1].id + 1;
		(*atoms)[currentIDAtoms].molType = 3;
		(*atoms)[currentIDAtoms].atomType = 10;
		(*atoms)[currentIDAtoms].charge = 0.435;
		(*atoms)[currentIDAtoms].x = butanediol[i].x + 3.67;
		(*atoms)[currentIDAtoms].y = butanediol[i].y - 1.80;
		(*atoms)[currentIDAtoms].z = butanediol[i].z + 0.92;
		currentIDAtoms++;

		// Adding bonds
		// O---H bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 8;
		(*bonds)[currentIDBonds].atom1 = H1_ID;
		(*bonds)[currentIDBonds].atom2 = O1_ID;
		currentIDBonds++;

		// O---C bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 9;
		(*bonds)[currentIDBonds].atom1 = O1_ID;
		(*bonds)[currentIDBonds].atom2 = C1_ID;
		currentIDBonds++;

		// C---C bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 10;
		(*bonds)[currentIDBonds].atom1 = C1_ID;
		(*bonds)[currentIDBonds].atom2 = C2_ID;
		currentIDBonds++;

		// C---C bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 10;
		(*bonds)[currentIDBonds].atom1 = C2_ID;
		(*bonds)[currentIDBonds].atom2 = C3_ID;
		currentIDBonds++;

		// C---C bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 10;
		(*bonds)[currentIDBonds].atom1 = C3_ID;
		(*bonds)[currentIDBonds].atom2 = C4_ID;
		currentIDBonds++;

		// C---O bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 9;
		(*bonds)[currentIDBonds].atom1 = C4_ID;
		(*bonds)[currentIDBonds].atom2 = O2_ID;
		currentIDBonds++;

		// O---H bond
		(*bonds)[currentIDBonds].id = (*bonds)[currentIDBonds-1].id + 1;
		(*bonds)[currentIDBonds].bondType = 8;
		(*bonds)[currentIDBonds].atom1 = O2_ID;
		(*bonds)[currentIDBonds].atom2 = H2_ID;
		currentIDBonds++;

		// Adding angles
		// H---O---C angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 11;
		(*angles)[currentIDAngles].atom1 = H1_ID;
		(*angles)[currentIDAngles].atom2 = O1_ID;
		(*angles)[currentIDAngles].atom3 = C1_ID;
		currentIDAngles++;

		// O---C---C angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 12;
		(*angles)[currentIDAngles].atom1 = O1_ID;
		(*angles)[currentIDAngles].atom2 = C1_ID;
		(*angles)[currentIDAngles].atom3 = C2_ID;
		currentIDAngles++;

		// C---C---C angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 13;
		(*angles)[currentIDAngles].atom1 = C1_ID;
		(*angles)[currentIDAngles].atom2 = C2_ID;
		(*angles)[currentIDAngles].atom3 = C3_ID;
		currentIDAngles++;

		// C---C---C angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 13;
		(*angles)[currentIDAngles].atom1 = C2_ID;
		(*angles)[currentIDAngles].atom2 = C3_ID;
		(*angles)[currentIDAngles].atom3 = C4_ID;
		currentIDAngles++;

		// C---C---O angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 12;
		(*angles)[currentIDAngles].atom1 = C3_ID;
		(*angles)[currentIDAngles].atom2 = C4_ID;
		(*angles)[currentIDAngles].atom3 = O2_ID;
		currentIDAngles++;

		// C---O---H angle
		(*angles)[currentIDAngles].id = (*angles)[currentIDAngles-1].id + 1;
		(*angles)[currentIDAngles].angleType = 11;
		(*angles)[currentIDAngles].atom1 = C4_ID;
		(*angles)[currentIDAngles].atom2 = O2_ID;
		(*angles)[currentIDAngles].atom3 = H2_ID;
		currentIDAngles++;

		// Adding dihedrals
		// H---O---C---C dihedral
		(*dihedrals)[currentIDDihedrals].id = (*dihedrals)[currentIDDihedrals-1].id + 1;
		(*dihedrals)[currentIDDihedrals].dihedralType = 10;
		(*dihedrals)[currentIDDihedrals].atom1 = H1_ID;
		(*dihedrals)[currentIDDihedrals].atom2 = O1_ID;
		(*dihedrals)[currentIDDihedrals].atom3 = C1_ID;
		(*dihedrals)[currentIDDihedrals].atom4 = C2_ID;
		currentIDDihedrals++;

		// O---C---C---C dihedral
		(*dihedrals)[currentIDDihedrals].id = (*dihedrals)[currentIDDihedrals-1].id + 1;
		(*dihedrals)[currentIDDihedrals].dihedralType = 11;
		(*dihedrals)[currentIDDihedrals].atom1 = O1_ID;
		(*dihedrals)[currentIDDihedrals].atom2 = C1_ID;
		(*dihedrals)[currentIDDihedrals].atom3 = C2_ID;
		(*dihedrals)[currentIDDihedrals].atom4 = C3_ID;
		currentIDDihedrals++;

		// C---C---C---C dihedral
		(*dihedrals)[currentIDDihedrals].id = (*dihedrals)[currentIDDihedrals-1].id + 1;
		(*dihedrals)[currentIDDihedrals].dihedralType = 12;
		(*dihedrals)[currentIDDihedrals].atom1 = C1_ID;
		(*dihedrals)[currentIDDihedrals].atom2 = C2_ID;
		(*dihedrals)[currentIDDihedrals].atom3 = C3_ID;
		(*dihedrals)[currentIDDihedrals].atom4 = C4_ID;
		currentIDDihedrals++;

		// C---C---C---O dihedral
		(*dihedrals)[currentIDDihedrals].id = (*dihedrals)[currentIDDihedrals-1].id + 1;
		(*dihedrals)[currentIDDihedrals].dihedralType = 11;
		(*dihedrals)[currentIDDihedrals].atom1 = C2_ID;
		(*dihedrals)[currentIDDihedrals].atom2 = C3_ID;
		(*dihedrals)[currentIDDihedrals].atom3 = C4_ID;
		(*dihedrals)[currentIDDihedrals].atom4 = O2_ID;
		currentIDDihedrals++;

		// C---C---O---H dihedral
		(*dihedrals)[currentIDDihedrals].id = (*dihedrals)[currentIDDihedrals-1].id + 1;
		(*dihedrals)[currentIDDihedrals].dihedralType = 10;
		(*dihedrals)[currentIDDihedrals].atom1 = C3_ID;
		(*dihedrals)[currentIDDihedrals].atom2 = C4_ID;
		(*dihedrals)[currentIDDihedrals].atom3 = O2_ID;
		(*dihedrals)[currentIDDihedrals].atom4 = H2_ID;
		currentIDDihedrals++;
	}
}

void recalculateNAtoms (DATA_ATOMS *atoms, DATAFILE_INFO *datafile)
{
	for (int i = 0; i < (*datafile).nAtoms; ++i)
	{
		if (atoms[i].atomType > (*datafile).nAtomTypes)
			(*datafile).nAtomTypes = atoms[i].atomType;
	}
}

void recalculateNBonds (DATA_BONDS *bonds, DATAFILE_INFO *datafile)
{
	for (int i = 0; i < (*datafile).nBonds; ++i)
	{
		if (bonds[i].bondType > (*datafile).nBondTypes)
			(*datafile).nBondTypes = bonds[i].bondType;
	}
}

void recalculateNAngles (DATA_ANGLES *angles, DATAFILE_INFO *datafile)
{
	for (int i = 0; i < (*datafile).nAngles; ++i)
	{
		if (angles[i].angleType > (*datafile).nAngleTypes)
			(*datafile).nAngleTypes = angles[i].angleType;
	}
}

void recalculateNDihedrals (DATA_DIHEDRALS *dihedrals, DATAFILE_INFO *datafile)
{
	for (int i = 0; i < (*datafile).nDihedrals; ++i)
	{
		if (dihedrals[i].dihedralType > (*datafile).nDihedralTypes)
			(*datafile).nDihedralTypes = dihedrals[i].dihedralType;
	}
}

void print_datafileHeader (FILE *output, BOUNDS simBoxDimension, DATAFILE_INFO datafile)
{
	fprintf(output, "Created by you v1.8.1 on today, this month, this year, current time.\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n%d impropers\n\n%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n%d improper types\n\n%.0f %.0f xlo xhi\n%.0f %.0f ylo yhi\n%.0f %.0f zlo zhi\n\nMasses\n\n1 1.008    # H of NaPSS\n2 12.011   # C of NaPSS\n3 13.018   # CH of NaPSS\n4 14.026   # CH2 of NaPSS\n5 32.065   # S of NaPSS\n6 15.999   # O of NaPSS\n7 22.989   # Na of NaPSS\n8 15.999   # O of H2O\n9 1.008    # H of H2O\n10 1.008   # H of butanediol\n11 15.999  # O of butanediol\n12 14.1707 # CH2_1 of butanediol\n\nAtoms\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers, datafile.nAtomTypes, datafile.nBondTypes, datafile.nAngleTypes, datafile.nDihedralTypes, datafile.nImproperTypes, simBoxDimension.xlo, simBoxDimension.xhi, simBoxDimension.ylo, simBoxDimension.yhi, simBoxDimension.zlo, simBoxDimension.zhi);
	fflush (output);
}

void print_dataAtoms (FILE *output, DATA_ATOMS *atoms, DATAFILE_INFO datafile)
{
	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		fprintf(output, "%d %d %d %.4f %.4f %.4f %.4f\n", atoms[i].id, atoms[i].molType, atoms[i].atomType, atoms[i].charge, atoms[i].x, atoms[i].y, atoms[i].z);
		fflush (output);
	}
}

void print_XYZatoms (FILE *outputXYZ, DATA_ATOMS *atoms, DATAFILE_INFO datafile)
{
	fprintf(outputXYZ, "%d\n", datafile.nAtoms);
	fprintf(outputXYZ, "%s\n", "#This XYZ file can be loaded after psf file for visualization.");

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		fprintf(outputXYZ, "C %.4f %.4f %.4f\n", atoms[i].x, atoms[i].y, atoms[i].z);
		fflush (outputXYZ);
	}
}

void print_dataBonds (FILE *output, DATA_BONDS *bonds, DATAFILE_INFO datafile)
{
	fprintf(output, "\nBonds\n\n");
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		fprintf(output, "%d %d %d %d\n", bonds[i].id, bonds[i].bondType, bonds[i].atom1, bonds[i].atom2);
		fflush (output);
	}
}

void print_dataAngles (FILE *output, DATA_ANGLES *angles, DATAFILE_INFO datafile)
{
	fprintf(output, "\nAngles\n\n");

	for (int i = 0; i < datafile.nAngles; ++i)
	{
		fprintf(output, "%d %d %d %d %d\n", angles[i].id, angles[i].angleType, angles[i].atom1, angles[i].atom2, angles[i].atom3);
		fflush (output);
	}
}

void print_dataDihedrals (FILE *output, DATA_DIHEDRALS *dihedrals, DATAFILE_INFO datafile)
{
	fprintf(output, "\nDihedrals\n\n");

	for (int i = 0; i < datafile.nDihedrals; ++i)
	{
		fprintf(output, "%d %d %d %d %d %d\n", dihedrals[i].id, dihedrals[i].dihedralType, dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4);
		fflush (output);
	}
}

void print_dataImpropers (FILE *output, DATA_IMPROPERS *impropers, DATAFILE_INFO datafile)
{
	fprintf(output, "\nImpropers\n\n");

	for (int i = 0; i < datafile.nImpropers; ++i)
	{
		fprintf(output, "%d %d %d %d %d %d\n", impropers[i].id, impropers[i].improperType, impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4);
		fflush (output);
	}

}

int main(int argc, char const *argv[])
{
	printf("\n");

	if (argc == 1)
	{
		printf("\nERROR: Insufficient arguments passed.\n\n   ARGS TO PASS:\n   ~~~~~~~~~~~~~\n\n * argv[0] = program\n * argv[1] = input dump file name\n * argv[2] = input data file\n\n");
		exit (1);
	}
	int nAtoms = getNatoms (argv[1]), lineCount = 0;

	char *pipeString, lineString[1000];
	pipeString = (char *) malloc (500 * sizeof (char));
	sprintf (pipeString, "tail -%d %s", nAtoms, argv[1]);

	FILE *input, *output, *input2, *outputXYZ;
	input = popen (pipeString, "r");
	input2 = fopen (argv[1], "r");
	output = fopen ("solvated.data", "w");
	outputXYZ = fopen ("solvated.xyz", "w");

	DUMP *traj, com, dimLow, dimHigh, boxLength, dumpLow, dumpHigh, chainDimension;
	BOUNDS dumpDimension, simBoxDimension;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	dumpDimension = getDumpDimension (input2);

	// float maxDimension = findMaxChainDimension (argv[1], nAtoms);
	chainDimension = findMaxChainDimension_XYZ (argv[1], nAtoms);

	traj = readLastFrame (pipeString, nAtoms, dumpDimension);

	com = findCOM (traj, nAtoms);
	printf("\ncalculating the center of mass from LAMMPS dump file...\n");
	printf(" com.x: %f\n com.y: %f\n com.z: %f\n", com.x, com.y, com.z);

	traj = recenterCoords (traj, com, dumpDimension, nAtoms);

	simBoxDimension = computeSimBoxDimension (chainDimension);

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafile, datafile_raw;

	datafile = readData (argv[2], &atoms, &bonds, &angles, &dihedrals, &impropers);

	// Replace the chain coords from the data file with the dump file
	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		atoms[i].x = traj[i].x;
		atoms[i].y = traj[i].y;
		atoms[i].z = traj[i].z;
	}

	// Create solvent molecules
	float solventDistance, butanediolFraction, waterFraction;
	DATA_ATOMS *water, water_com, recenterDistance, *butanediol, butanediol_com;
	water_com.x = 0; water_com.y = 0; water_com.z = 0; butanediol_com.x = 0; butanediol_com.y = 0; butanediol_com.z = 0;

	printf("Enter the volume fraction of butanediol [0 to 1]: "); scanf ("%f", &butanediolFraction); printf("\n");
	int nButanediol = findNButanediol (butanediolFraction, simBoxDimension);
	waterFraction = 1.0 - butanediolFraction;
	int nWater = findNWater (waterFraction, simBoxDimension);

	water = (DATA_ATOMS *) calloc (nWater, sizeof (DATA_ATOMS));
	butanediol = (DATA_ATOMS *) calloc (nButanediol, sizeof (DATA_ATOMS));

	butanediol = populateButanediol (simBoxDimension, nButanediol);
	water = populateWater (simBoxDimension, nWater);

	// Creating topology information for water and butanediol
	// Reallocating memory for extra solvent atoms
	datafile_raw = datafile;
	datafile.nAtoms += (nWater * 3) + (nButanediol * 8);
	atoms = (DATA_ATOMS *) realloc (atoms, datafile.nAtoms * sizeof (DATA_ATOMS));

	datafile.nBonds += (nWater * 2) + (nButanediol * 7);
	bonds = (DATA_BONDS *) realloc (bonds, datafile.nBonds * sizeof (DATA_BONDS));

	datafile.nAngles += nWater + (nButanediol * 6);
	angles = (DATA_ANGLES *) realloc (angles, datafile.nAngles * sizeof (DATA_ANGLES));

	datafile.nDihedrals += (nButanediol * 5);
	dihedrals = (DATA_DIHEDRALS *) realloc (dihedrals, datafile.nDihedrals * sizeof (DATA_DIHEDRALS));

	createButanediolTopology (&atoms, &bonds, &angles, &dihedrals, nButanediol, butanediol, datafile_raw, datafile);
	createWaterTopology (&atoms, &bonds, &angles, nWater, water, datafile_raw, datafile);

	// Recalculating nAtomTypes, nBondsTypes, nAngleTypes and nDihedralTypes after adding new solvent molecules
	recalculateNAtoms (atoms, &datafile);
	recalculateNBonds (bonds, &datafile);
	recalculateNAngles (angles, &datafile);
	recalculateNDihedrals (dihedrals, &datafile);

	// Printing datafile
	print_datafileHeader (output, simBoxDimension, datafile);
	print_dataAtoms (output, atoms, datafile);
	print_XYZatoms (outputXYZ, atoms, datafile);
	print_dataBonds (output, bonds, datafile);
	print_dataAngles (output, angles, datafile);
	print_dataDihedrals (output, dihedrals, datafile);
	print_dataImpropers (output, impropers, datafile);


	// Zero the water and butanediol coordinates
	// for (int i = 0; i < nWater; ++i)
	// {
	// 	water[i].x = 0; water[i].y = 0; water[i].z = 0;
	// }

	// for (int i = 0; i < nButanediol; ++i)
	// {
	// 	butanediol[i].x = 0; butanediol[i].y = 0; butanediol[i].z = 0;
	// }

	// int nBins_x_water = (int) boxLength.x / 3.103, nBins_y_water = (int) boxLength.y / 3.103, nBins_z_water = (int) boxLength.z / 3.103, nBins_x_butanediol = (int) boxLength.x / 8.0, nBins_y_butanediol = (int) boxLength.y / 5.0, nBins_z_butanediol = (int) boxLength.z / 3.0;

	// if (nWater == 0)
	// {
	// 	nBins_x_water = 0;
	// 	nBins_y_water = 0;
	// 	nBins_z_water = 0;
	// }

	// if (nButanediol == 0)
	// {
	// 	nBins_x_butanediol = 0;
	// 	nBins_y_butanediol = 0;
	// 	nBins_z_butanediol = 0;
	// }

	// printf("nBins_x_water: %d; nBins_y_water: %d; nBins_z_water: %d\n", nBins_x_water, nBins_y_water, nBins_z_water);
	// printf("nBins_x_butanediol: %d; nBins_y_butanediol: %d; nBins_z_butanediol: %d\n", nBins_x_butanediol, nBins_y_butanediol, nBins_z_butanediol);
	// printf("Total number of water molecules that will fit in the simulation box: %d\n", nBins_x_water * nBins_y_water * nBins_z_water);

	// int nWater_remaining = nWater, nButanediol_remaining = nButanediol, isRandomSolvent = atoi (argv[5]);
	// int currentWater = 0, currentButanediol = 0;
	// int butanediolOverlap = 1, waterOverlap = 1;

	// To Do: 
	// Don't ask for input on number of water molecules to be added.
	// nButanediol can be user defined input
	// nWater is based on the remaining free volume
	// Convert these into functions and re-write them properly
	// In long term, write a generalized code in python

	// Butanediol molecules are added first.
	// Butanediol molecules are added only after checking the distance
	// between the new molecules and the existing NaPSS chain
	// for (int i = 0; i < nBins_x_butanediol; ++i)
	// {
	// 	for (int j = 0; j < nBins_y_butanediol; ++j)
	// 	{
	// 		for (int k = 0; k < nBins_z_butanediol; ++k)
	// 		{
	// 			if (nButanediol_remaining)
	// 			{
	// 				butanediol[currentButanediol].x = ((i + 1) * 8) - (8 / 2);
	// 				butanediol[currentButanediol].y = ((j + 1) * 5) - (5 / 2);
	// 				butanediol[currentButanediol].z = ((k + 1) * 3) - (3 / 2);

	// 				// Checking the distance between newly added butanediol and the polymer chain
	// 				for (int j = 0; j < lineCount; ++j)
	// 				{
	// 					solventDistance = sqrt (pow ((butanediol[currentButanediol].x - traj[j].x), 2) + pow ((butanediol[currentButanediol].y - traj[j].y), 2) + pow ((butanediol[currentButanediol].z - traj[j].z), 2));

	// 					if (solventDistance < 9)
	// 						butanediolOverlap = 1;
	// 					else
	// 						butanediolOverlap = 0;
	// 				}

	// 				// Indices are updated only if there is no overlap
	// 				if (butanediolOverlap == 0)
	// 				{
	// 					currentButanediol++;
	// 					nButanediol_remaining--;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	// // Water molecules are added after butanediol
	// // Water molecules are only added after checking the distance between 
	// // the new molecules and the existing NaPSS chain, and butanediol molecules
	// for (int i = 0; i < nBins_x_water; ++i)
	// {
	// 	for (int j = 0; j < nBins_y_water; ++j)
	// 	{
	// 		for (int k = 0; k < nBins_z_water; ++k)
	// 		{
	// 			if (nWater_remaining)
	// 			{
	// 				water[currentWater].x = ((i + 1) * 3.103) - (3.103 / 2);
	// 				water[currentWater].y = ((j + 1) * 3.103) - (3.103 / 2);
	// 				water[currentWater].z = ((k + 1) * 3.103) - (3.103 / 2);

	// 				// Checking the distance between newly added water and the polymer chain
	// 				for (int j = 0; j < lineCount; ++j)
	// 				{
	// 					solventDistance = sqrt (pow ((water[currentWater].x - traj[j].x), 2) + pow ((water[currentWater].y - traj[j].y), 2) + pow ((water[currentWater].z - traj[j].z), 2));

	// 					if (solventDistance < 3.5)
	// 						waterOverlap = 1;
	// 					else
	// 						waterOverlap = 0;
	// 				}
	// 				// Checking the distance between newly added water and butanediol molecules
	// 				for (int j = 0; j < nButanediol; ++j)
	// 				{
	// 					if (butanediol[i].x > 0 || butanediol[i].y > 0 || butanediol[i].z > 0)
	// 					{
	// 						solventDistance = sqrt (pow ((water[currentWater].x - butanediol[j].x), 2) + pow ((water[currentWater].y - butanediol[j].y), 2) + pow ((water[currentWater].z - butanediol[j].z), 2));

	// 						if (solventDistance < 9)
	// 							waterOverlap = 1;
	// 					}
	// 				}

	// 				// Indices are updated only if there is no overlap
	// 				if (waterOverlap == 0)
	// 				{
	// 					currentWater++;
	// 					nWater_remaining--;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	// nWater = currentWater;
	// nButanediol = currentButanediol;

	// printf("Total number of water molecules added: %d\n", nWater);
	// printf("Total number of butanediol molecules added: %d\n", nButanediol);

	// Finding the center of mass
	// for (int i = 0; i < nWater; ++i)
	// {
	// 	water_com.x += water[i].x; water_com.y += water[i].y; water_com.z += water[i].z;
	// }

	// for (int i = 0; i < nButanediol; ++i)
	// {
	// 	butanediol_com.x += butanediol[i].x; butanediol_com.y += butanediol[i].y; butanediol_com.z += butanediol[i].z;
	// }

	// water_com.x /= nWater; water_com.y /= nWater; water_com.z /= nWater; butanediol_com.x /= nButanediol; butanediol_com.y /= nButanediol; butanediol_com.z /= nButanediol;

	// printf("Calculating center of mass for water molecules\n water_com.x: %f; water_com.y: %f; water_com.z: %f\n", water_com.x, water_com.y, water_com.z);
	// printf("Calculating center of mass for butanediol molecules\n butanediol_com.x: %f; butanediol_com.y: %f; butanediol_com.z: %f\n", butanediol_com.x, butanediol_com.y, butanediol_com.z);

	// Recenter solvent molecules
	// for (int i = 0; i < nWater; ++i)
	// {
	// 	water[i].x -= water_com.x; water[i].y -= water_com.y; water[i].z -= water_com.z;
	// }

	// for (int i = 0; i < nButanediol; ++i)
	// {
	// 	butanediol[i].x -= butanediol_com.x; butanediol[i].y -= butanediol_com.y; butanediol[i].z -= butanediol_com.z;
	// }

	// Recalculating the center of mass for all solvent molecules
	// water_com.x = 0; water_com.y = 0; water_com.z = 0;
	// for (int i = 0; i < nWater; ++i)
	// {
	// 	water_com.x += water[i].x; water_com.y += water[i].y; water_com.z += water[i].z;
	// }

	// for (int i = 0; i < nButanediol; ++i)
	// {
	// 	butanediol_com.x += butanediol[i].x; butanediol_com.y += butanediol[i].y; butanediol_com.z += butanediol[i].z;
	// }

	// water_com.x /= nWater; water_com.y /= nWater; water_com.z /= nWater;
	// butanediol_com.x /= nButanediol; butanediol_com.y /= nButanediol; butanediol_com.z /= nButanediol;

	// printf("\nRecentering water molecules\n water_com.x: %f; water_com.y: %f; water_com.z: %f\n", water_com.x, water_com.y, water_com.z);
	// printf("\nRecentering butanediol molecules\n butanediol_com.x: %f; butanediol_com.y: %f; butanediol_com.z: %f\n", butanediol_com.x, butanediol_com.y, butanediol_com.z);

	// Create topology for water molecules
	// int newNAtoms = datafile.nAtoms + (nWater * 3) + (nButanediol * 8), currentIDAtoms = datafile.nAtoms;
	// atoms = (DATA_ATOMS *) realloc (atoms, newNAtoms * sizeof (DATA_ATOMS));

	// int newNBonds = datafile.nBonds + (nWater * 2) + (nButanediol * 7), currentIDBonds = datafile.nBonds;
	// bonds = (DATA_BONDS *) realloc (bonds, newNBonds * sizeof (DATA_BONDS));

	// int newNAngles = datafile.nAngles + nWater + (nButanediol * 6), currentIDAngles = datafile.nAngles;
	// angles = (DATA_ANGLES *) realloc (angles, newNAngles * sizeof (DATA_ANGLES));

	// int newNDihedrals = datafile.nDihedrals + (nButanediol * 5), currentIDDihedrals = datafile.nDihedrals;
	// dihedrals = (DATA_DIHEDRALS *) realloc (dihedrals, newNDihedrals * sizeof (DATA_DIHEDRALS));

	// for (int i = 0; i < nWater; ++i)
	// {
	// 	// Adding oxygen molecule
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 2;
	// 	atoms[currentIDAtoms].atomType = 8;
	// 	atoms[currentIDAtoms].charge = -0.83;
	// 	atoms[currentIDAtoms].x = water[i].x;
	// 	atoms[currentIDAtoms].y = water[i].y;
	// 	atoms[currentIDAtoms].z = water[i].z;
	// 	currentIDAtoms++;

	// 	// Adding hydrogens
	// 	// Hydrogen 1
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 2;
	// 	atoms[currentIDAtoms].atomType = 9;
	// 	atoms[currentIDAtoms].charge = 0.415;
	// 	atoms[currentIDAtoms].x = water[i].x - 0.32;
	// 	atoms[currentIDAtoms].y = water[i].y + 0.13;
	// 	atoms[currentIDAtoms].z = water[i].z - 0.91;
	// 	currentIDAtoms++;

	// 	// Hydrogen 2
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 2;
	// 	atoms[currentIDAtoms].atomType = 9;
	// 	atoms[currentIDAtoms].charge = 0.415;
	// 	atoms[currentIDAtoms].x = water[i].x + 0.97;
	// 	atoms[currentIDAtoms].y = water[i].y;
	// 	atoms[currentIDAtoms].z = water[i].z;
	// 	currentIDAtoms++;

	// 	// Adding bonds (bond1 and bond2 are identical)
	// 	// Bond 1
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 7;
	// 	bonds[currentIDBonds].atom1 = atoms[currentIDAtoms-1].id;
	// 	bonds[currentIDBonds].atom2 = atoms[currentIDAtoms-3].id;
	// 	currentIDBonds++;

	// 	// Bond 2
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 7;
	// 	bonds[currentIDBonds].atom1 = atoms[currentIDAtoms-2].id;
	// 	bonds[currentIDBonds].atom2 = atoms[currentIDAtoms-3].id;
	// 	currentIDBonds++;

	// 	// Adding angles
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 10;
	// 	angles[currentIDAngles].atom1 = atoms[currentIDAtoms-1].id;
	// 	angles[currentIDAngles].atom2 = atoms[currentIDAtoms-3].id;
	// 	angles[currentIDAngles].atom3 = atoms[currentIDAtoms-2].id;
	// 	currentIDAngles++;
	// }

	// int H1id, O1id, C1id, C2id, C3id, C4id, O2id, H2id;

	// for (int i = 0; i < nButanediol; ++i)
	// {
	// 	// Defined by TraPPE force field
	// 	// Adding atomic coordinates
	// 	// H
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	H1id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 10;
	// 	atoms[currentIDAtoms].charge = 0.435;
	// 	atoms[currentIDAtoms].x = butanediol[i].x - 2.18;
	// 	atoms[currentIDAtoms].y = butanediol[i].y + 1.29;
	// 	atoms[currentIDAtoms].z = butanediol[i].z - 1.26;
	// 	currentIDAtoms++;

	// 	// O
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	O1id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 11;
	// 	atoms[currentIDAtoms].charge = -0.7;
	// 	atoms[currentIDAtoms].x = butanediol[i].x - 1.93;
	// 	atoms[currentIDAtoms].y = butanediol[i].y + 0.34;
	// 	atoms[currentIDAtoms].z = butanediol[i].z - 1.36;
	// 	currentIDAtoms++;

	// 	// CH2
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	C1id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 12;
	// 	atoms[currentIDAtoms].charge = 0.265;
	// 	atoms[currentIDAtoms].x = butanediol[i].x - 0.53;
	// 	atoms[currentIDAtoms].y = butanediol[i].y + 0.34;
	// 	atoms[currentIDAtoms].z = butanediol[i].z - 1.37;
	// 	currentIDAtoms++;

	// 	// CH2
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	C2id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 12;
	// 	atoms[currentIDAtoms].charge = 0.0;
	// 	atoms[currentIDAtoms].x = butanediol[i].x;
	// 	atoms[currentIDAtoms].y = butanediol[i].y;
	// 	atoms[currentIDAtoms].z = butanediol[i].z;
	// 	currentIDAtoms++;

	// 	// CH2
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	C3id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 12;
	// 	atoms[currentIDAtoms].charge = 0.0;
	// 	atoms[currentIDAtoms].x = butanediol[i].x + 1.52;
	// 	atoms[currentIDAtoms].y = butanediol[i].y;
	// 	atoms[currentIDAtoms].z = butanediol[i].z;
	// 	currentIDAtoms++;

	// 	// CH2
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	C4id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 12;
	// 	atoms[currentIDAtoms].charge = 0.265;
	// 	atoms[currentIDAtoms].x = butanediol[i].x + 2.05;
	// 	atoms[currentIDAtoms].y = butanediol[i].y - 1.25;
	// 	atoms[currentIDAtoms].z = butanediol[i].z - 0.68;
	// 	currentIDAtoms++;

	// 	// O
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	O2id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 11;
	// 	atoms[currentIDAtoms].charge = -0.7;
	// 	atoms[currentIDAtoms].x = butanediol[i].x + 3.44;
	// 	atoms[currentIDAtoms].y = butanediol[i].y - 1.23;
	// 	atoms[currentIDAtoms].z = butanediol[i].z - 0.67;
	// 	currentIDAtoms++;

	// 	// H
	// 	atoms[currentIDAtoms].id = atoms[currentIDAtoms-1].id + 1;
	// 	H2id = atoms[currentIDAtoms-1].id + 1;
	// 	atoms[currentIDAtoms].molType = 3;
	// 	atoms[currentIDAtoms].atomType = 10;
	// 	atoms[currentIDAtoms].charge = 0.435;
	// 	atoms[currentIDAtoms].x = butanediol[i].x + 3.67;
	// 	atoms[currentIDAtoms].y = butanediol[i].y - 1.80;
	// 	atoms[currentIDAtoms].z = butanediol[i].z + 0.92;
	// 	currentIDAtoms++;

	// 	// Adding bonds
	// 	// O---H bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 8;
	// 	bonds[currentIDBonds].atom1 = H1id;
	// 	bonds[currentIDBonds].atom2 = O1id;
	// 	currentIDBonds++;

	// 	// O---C bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 9;
	// 	bonds[currentIDBonds].atom1 = O1id;
	// 	bonds[currentIDBonds].atom2 = C1id;
	// 	currentIDBonds++;

	// 	// C---C bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 10;
	// 	bonds[currentIDBonds].atom1 = C1id;
	// 	bonds[currentIDBonds].atom2 = C2id;
	// 	currentIDBonds++;

	// 	// C---C bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 10;
	// 	bonds[currentIDBonds].atom1 = C2id;
	// 	bonds[currentIDBonds].atom2 = C3id;
	// 	currentIDBonds++;

	// 	// C---C bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 10;
	// 	bonds[currentIDBonds].atom1 = C3id;
	// 	bonds[currentIDBonds].atom2 = C4id;
	// 	currentIDBonds++;

	// 	// C---O bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 9;
	// 	bonds[currentIDBonds].atom1 = C4id;
	// 	bonds[currentIDBonds].atom2 = O2id;
	// 	currentIDBonds++;

	// 	// O---H bond
	// 	bonds[currentIDBonds].id = bonds[currentIDBonds-1].id + 1;
	// 	bonds[currentIDBonds].bondType = 8;
	// 	bonds[currentIDBonds].atom1 = O2id;
	// 	bonds[currentIDBonds].atom2 = H2id;
	// 	currentIDBonds++;

	// 	// Adding angles
	// 	// H---O---C angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 11;
	// 	angles[currentIDAngles].atom1 = H1id;
	// 	angles[currentIDAngles].atom2 = O1id;
	// 	angles[currentIDAngles].atom3 = C1id;
	// 	currentIDAngles++;

	// 	// O---C---C angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 12;
	// 	angles[currentIDAngles].atom1 = O1id;
	// 	angles[currentIDAngles].atom2 = C1id;
	// 	angles[currentIDAngles].atom3 = C2id;
	// 	currentIDAngles++;

	// 	// C---C---C angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 13;
	// 	angles[currentIDAngles].atom1 = C1id;
	// 	angles[currentIDAngles].atom2 = C2id;
	// 	angles[currentIDAngles].atom3 = C3id;
	// 	currentIDAngles++;

	// 	// C---C---C angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 13;
	// 	angles[currentIDAngles].atom1 = C2id;
	// 	angles[currentIDAngles].atom2 = C3id;
	// 	angles[currentIDAngles].atom3 = C4id;
	// 	currentIDAngles++;

	// 	// C---C---O angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 12;
	// 	angles[currentIDAngles].atom1 = C3id;
	// 	angles[currentIDAngles].atom2 = C4id;
	// 	angles[currentIDAngles].atom3 = O2id;
	// 	currentIDAngles++;

	// 	// C---O---H angle
	// 	angles[currentIDAngles].id = angles[currentIDAngles-1].id + 1;
	// 	angles[currentIDAngles].angleType = 11;
	// 	angles[currentIDAngles].atom1 = C4id;
	// 	angles[currentIDAngles].atom2 = O2id;
	// 	angles[currentIDAngles].atom3 = H2id;
	// 	currentIDAngles++;

	// 	// Adding dihedrals
	// 	// H---O---C---C dihedral
	// 	dihedrals[currentIDDihedrals].id = dihedrals[currentIDDihedrals-1].id + 1;
	// 	dihedrals[currentIDDihedrals].dihedralType = 10;
	// 	dihedrals[currentIDDihedrals].atom1 = H1id;
	// 	dihedrals[currentIDDihedrals].atom2 = O1id;
	// 	dihedrals[currentIDDihedrals].atom3 = C1id;
	// 	dihedrals[currentIDDihedrals].atom4 = C2id;
	// 	currentIDDihedrals++;

	// 	// O---C---C---C dihedral
	// 	dihedrals[currentIDDihedrals].id = dihedrals[currentIDDihedrals-1].id + 1;
	// 	dihedrals[currentIDDihedrals].dihedralType = 11;
	// 	dihedrals[currentIDDihedrals].atom1 = O1id;
	// 	dihedrals[currentIDDihedrals].atom2 = C1id;
	// 	dihedrals[currentIDDihedrals].atom3 = C2id;
	// 	dihedrals[currentIDDihedrals].atom4 = C3id;
	// 	currentIDDihedrals++;

	// 	// C---C---C---C dihedral
	// 	dihedrals[currentIDDihedrals].id = dihedrals[currentIDDihedrals-1].id + 1;
	// 	dihedrals[currentIDDihedrals].dihedralType = 12;
	// 	dihedrals[currentIDDihedrals].atom1 = C1id;
	// 	dihedrals[currentIDDihedrals].atom2 = C2id;
	// 	dihedrals[currentIDDihedrals].atom3 = C3id;
	// 	dihedrals[currentIDDihedrals].atom4 = C4id;
	// 	currentIDDihedrals++;

	// 	// C---C---C---O dihedral
	// 	dihedrals[currentIDDihedrals].id = dihedrals[currentIDDihedrals-1].id + 1;
	// 	dihedrals[currentIDDihedrals].dihedralType = 11;
	// 	dihedrals[currentIDDihedrals].atom1 = C2id;
	// 	dihedrals[currentIDDihedrals].atom2 = C3id;
	// 	dihedrals[currentIDDihedrals].atom3 = C4id;
	// 	dihedrals[currentIDDihedrals].atom4 = O2id;
	// 	currentIDDihedrals++;

	// 	// C---C---O---H dihedral
	// 	dihedrals[currentIDDihedrals].id = dihedrals[currentIDDihedrals-1].id + 1;
	// 	dihedrals[currentIDDihedrals].dihedralType = 10;
	// 	dihedrals[currentIDDihedrals].atom1 = C3id;
	// 	dihedrals[currentIDDihedrals].atom2 = C4id;
	// 	dihedrals[currentIDDihedrals].atom3 = O2id;
	// 	dihedrals[currentIDDihedrals].atom4 = H2id;
	// 	currentIDDihedrals++;
	// }

	// Recalculting the number of atom, bond, angle, dihedral and improper types
	// for (int i = 0; i < newNAtoms; ++i)
	// {
	// 	if (atoms[i].atomType > datafile.nAtomTypes)
	// 		datafile.nAtomTypes = atoms[i].atomType;
	// }

	// for (int i = 0; i < newNBonds; ++i)
	// {
	// 	if (bonds[i].bondType > datafile.nBondTypes)
	// 		datafile.nBondTypes = bonds[i].bondType;
	// }

	// for (int i = 0; i < newNAngles; ++i)
	// {
	// 	if (angles[i].angleType > datafile.nAngleTypes)
	// 		datafile.nAngleTypes = angles[i].angleType;
	// }

	// for (int i = 0; i < newNDihedrals; ++i)
	// {
	// 	if (dihedrals[i].dihedralType > datafile.nDihedralTypes)
	// 		datafile.nDihedralTypes = dihedrals[i].dihedralType;
	// }

	// for (int i = 0; i < datafile.nImpropers; ++i)
	// {
	// 	if (impropers[i].improperType > datafile.nImproperTypes)
	// 		datafile.nImproperTypes = impropers[i].improperType;
	// }

	// printf("\nRecalculated values for\n\n Atom types: %d\n Bond types: %d\n Angle types: %d\n Dihedral types: %d\n Improper types: %d\n", datafile.nAtomTypes, datafile.nBondTypes, datafile.nAngleTypes, datafile.nDihedralTypes, datafile.nImproperTypes);

	// Print the final data file
	// fprintf(output, "Created by you v1.8.1 on today, this month, this year, current time.\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n%d impropers\n\n%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n%d improper types\n\n%.0f %.0f xlo xhi\n%.0f %.0f ylo yhi\n%.0f %.0f zlo zhi\n\nMasses\n\n1 1.008    # H of NaPSS\n2 12.011   # C of NaPSS\n3 13.018   # CH of NaPSS\n4 14.026   # CH2 of NaPSS\n5 32.065   # S of NaPSS\n6 15.999   # O of NaPSS\n7 22.989   # Na of NaPSS\n8 15.999   # O of H2O\n9 1.008    # H of H2O\n10 1.008   # H of butanediol\n11 15.999  # O of butanediol\n12 14.1707 # CH2_1 of butanediol\n\nAtoms\n\n", newNAtoms, newNBonds, newNAngles, newNDihedrals, datafile.nImpropers, datafile.nAtomTypes, datafile.nBondTypes, datafile.nAngleTypes, datafile.nDihedralTypes, datafile.nImproperTypes, dimLow.x - 1, dimHigh.x + 1, dimLow.y - 1, dimHigh.y + 1, dimLow.z - 1, dimHigh.z + 1);

	// 1 1.008    # H
	// 2 12.011   # C
	// 3 13.018   # CH
	// 4 14.026   # CH2
	// 5 32.065   # S
	// 6 15.999   # O
	// 7 22.989   # Na
	// 8 15.999   # O from H2O
	// 9 1.008    # H from H2O
	// 10 1.008   # H from butanediol
	// 11 15.999  # O from butanediol
	// 12 14.1707 # CH2_1 from butanediol
	// 13 14.1707 # CH2_2 from butanediol

	// fprintf(outputXYZ, "%d\n", newNAtoms);
	// fprintf(outputXYZ, "%s\n", "#This XYZ file can be loaded after psf file for visualization.");
	// for (int i = 0; i < newNAtoms; ++i)
	// {
	// 	fprintf(output, "%d %d %d %.4f %.4f %.4f %.4f\n", atoms[i].id, atoms[i].molType, atoms[i].atomType, atoms[i].charge, atoms[i].x, atoms[i].y, atoms[i].z);
	// 	fprintf(outputXYZ, "C %.4f %.4f %.4f\n", atoms[i].x, atoms[i].y, atoms[i].z);
	// }

	// fprintf(output, "\nBonds\n\n");

	// for (int i = 0; i < newNBonds; ++i)
	// 	fprintf(output, "%d %d %d %d\n", bonds[i].id, bonds[i].bondType, bonds[i].atom1, bonds[i].atom2);

	// fprintf(output, "\nAngles\n\n");

	// for (int i = 0; i < newNAngles; ++i)
	// 	fprintf(output, "%d %d %d %d %d\n", angles[i].id, angles[i].angleType, angles[i].atom1, angles[i].atom2, angles[i].atom3);

	// fprintf(output, "\nDihedrals\n\n");

	// for (int i = 0; i < newNDihedrals; ++i)
	// 	fprintf(output, "%d %d %d %d %d %d\n", dihedrals[i].id, dihedrals[i].dihedralType, dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4);

	// fprintf(output, "\nImpropers\n\n");

	// for (int i = 0; i < datafile.nImpropers; ++i)
	// 	fprintf(output, "%d %d %d %d %d %d\n", impropers[i].id, impropers[i].improperType, impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4);

	fclose (input);
	fclose (output);
	fclose (input2);
	fclose (outputXYZ);

	return 0;
}