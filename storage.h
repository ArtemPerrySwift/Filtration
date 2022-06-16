#pragma once
#include<string>
#include<fstream>

struct Storage
{
//	StorageData* storageData;
	
	int count;

	/// <summary>
	/// ���������� ������ �� ���������� �����
	/// </summary>
	/// <param name="path">���� � �����</param>
	/// <param name="count">���������� ���������, ������� ����� �������</param>
	/// <returns></returns>
	virtual bool read(std::string path, int count = -1);

	/// <summary>
	/// ���������� ������ �� ��������� �����
	/// </summary>
	/// <param name="path"></param>
	/// <param name="count"></param>
	/// <returns></returns>
	virtual bool readBin(std::string path, int count = -1);

	virtual bool read(std::ifstream &in, std::string path, int count = -1);

protected:

	/// <summary>
	/// ��������� ������ ��� ���������
	/// </summary>
	/// <returns></returns>
	virtual bool allocateMemory() = 0;

	/// <summary>
	/// ���������� i-�� �������� ���������
	/// </summary>
	/// <param name="in"></param>
	/// <param name="i"></param>
	/// <returns></returns>
	virtual bool readData(std::ifstream& in, int i) = 0;
};


struct StorageData
{
	virtual bool readData(std::ifstream& in) = 0;
};

