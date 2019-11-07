#ifndef __GX_MATRIX_H
#define __GX_MATRIX_H

#include <vector>
#include <string>

template <typename T> class GXMatrix {
 private:
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

 public:
  GXMatrix();
  GXMatrix(unsigned _rows, unsigned _cols, const T& _initial);
  GXMatrix(const GXMatrix<T>& rhs);
  ~GXMatrix();

   // Increase size
  void add_row(std::vector<T> row);
  void add_row(unsigned n_rows, const T& val);

  // Operator overloading, for "standard" mathematical matrix operations                                                                                                                                                          
  GXMatrix<T>& operator=(const GXMatrix<T>& rhs);

  // Matrix mathematical operations                                                                                                                                                                                               
  GXMatrix<T> operator+(const GXMatrix<T>& rhs);
  GXMatrix<T>& operator+=(const GXMatrix<T>& rhs);
  GXMatrix<T> operator-(const GXMatrix<T>& rhs);
  GXMatrix<T>& operator-=(const GXMatrix<T>& rhs);
  GXMatrix<T> operator*(const GXMatrix<T>& rhs);
  GXMatrix<T>& operator*=(const GXMatrix<T>& rhs);
  GXMatrix<T> operator/(const GXMatrix<T>& rhs);
  GXMatrix<T>& operator/=(const GXMatrix<T>& rhs);
  
  GXMatrix<T> transpose();

  // Matrix/scalar operations                                                                                                                                                                                                     
  GXMatrix<T> operator+(const T& rhs);
  GXMatrix<T> operator-(const T& rhs);
  GXMatrix<T> operator*(const T& rhs);
  GXMatrix<T> operator/(const T& rhs);

  // Matrix/vector operations                                                                                                                                                                                                     
  //std::vector<T> operator*(const std::vector<T>& rhs);
  //std::vector<T> diag_vec();

  // Access the individual elements 
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;
  std::vector<T>& getRow(const unsigned& row);
  const std::vector<T>& getRow(const unsigned& row) const;
  // Access the row and column sizes
  unsigned get_rows() const;
  unsigned get_cols() const;
  std::vector<std::vector<T> > getMat();
  std::string toString();
};





//Formerly in .cpp
template <typename T>  GXMatrix<T>::GXMatrix(){
	this->rows = 0;
	this->cols = 0;
}


template <typename T>  GXMatrix<T>::GXMatrix(unsigned _rows, unsigned _cols, const T& _initial){
	this->rows = _rows;
	this->cols = _cols;
	for(int i = 0; i< _rows; i++){
		this->mat.push_back(std::vector<T>(_cols, _initial));
	}
}

template <typename T> GXMatrix<T>::GXMatrix(const GXMatrix<T>& rhs){
	this->rows = rhs.get_rows();
	this->cols = rhs.get_cols();
	std::vector<T> aux;
	for(unsigned i = 0; i < rows; i++){
		aux = std::vector<T>(rhs.getRow(i));
		this->mat.push_back(aux);
	}
}


template <typename T> GXMatrix<T>::~GXMatrix(){
}



template <typename T>  void GXMatrix<T>::add_row(std::vector<T> row){
	if(row.size() == rows){
		this->mat.push_back(row);
		this->rows++;
	}
}

template <typename T>  void GXMatrix<T>::add_row(unsigned n_rows, const T& _initial){
	for(int i = 0; i < n_rows; i++){
		this->mat.push_back(std::vector<T>(this->cols, _initial));
	}
}

  // Operator overloading, for "standard" mathematical matrix operations
                                                                                                                                                          
template <typename T>  GXMatrix<T>& GXMatrix<T>::operator=(const GXMatrix<T>& rhs){
	this->mat = rhs.mat;
	this->cols = rhs.get_cols();
	this->rows = rhs.get_rows();
	return *this;
}

  // Matrix mathematical operations
template <typename T>  GXMatrix<T> GXMatrix<T>::operator+(const GXMatrix<T>& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) += rhs(i,j);
	return result;
}
template <typename T>  GXMatrix<T>& GXMatrix<T>::operator+=(const GXMatrix<T>& rhs){
	*this = *this + rhs;
	return *this;
}
template <typename T>  GXMatrix<T> GXMatrix<T>::operator-(const GXMatrix<T>& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) -= rhs(i,j);
	return result;
}
template <typename T> GXMatrix<T>& GXMatrix<T>::operator-=(const GXMatrix<T>& rhs){
	*this = *this - rhs;
	return *this;
}
template <typename T>  GXMatrix<T> GXMatrix<T>::operator*(const GXMatrix<T>& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) *= rhs(i,j);
	return result;
}
template <typename T>  GXMatrix<T>& GXMatrix<T>::operator*=(const GXMatrix<T>& rhs){
	*this = *this * rhs;
	return *this;
}
template <typename T>  GXMatrix<T> GXMatrix<T>::operator/(const GXMatrix<T>& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) /= rhs(i,j);
	return result;
}
template <typename T>  GXMatrix<T>& GXMatrix<T>::operator/=(const GXMatrix<T>& rhs){
	*this = *this / rhs;
	return *this;
}

template <typename T> GXMatrix<T> GXMatrix<T>::transpose(){
	GXMatrix<T> aux(*this);
	GXMatrix<T> result(this->get_cols(), this->get_rows(), aux(0,0));
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) = aux(j,i);
	return result;
}

  // Matrix/scalar operations                                                                                                                                                                                                     
template <typename T> GXMatrix<T> GXMatrix<T>::operator+(const T& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) += rhs;
	return result;
}
template <typename T>  GXMatrix<T> GXMatrix<T>::operator-(const T& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) -= rhs;
	return result;
}
template <typename T> GXMatrix<T> GXMatrix<T>::operator*(const T& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) *= rhs;
	return result;
}
template <typename T> GXMatrix<T> GXMatrix<T>::operator/(const T& rhs){
	GXMatrix<T> result(*this);
	for(int i = 0; i < result.get_rows(); i++) for(int j = 0; j < result.get_cols(); j++) result(i, j) /= rhs;
	return result;
}

  // Matrix/vector operations                                                                                                                                                                                                     
  //std::vector<T> operator*(const std::vector<T>& rhs);
  //std::vector<T> diag_vec();

  // Access the individual elements      
                                                                                                                                                                                         
template <typename T> T& GXMatrix<T>::operator()(const unsigned& row, const unsigned& col){
	return this->mat[row][col];
}

template <typename T> const T& GXMatrix<T>::operator()(const unsigned& row, const unsigned& col) const{
	return this->mat[row][col];
}

template <typename T> const std::vector<T>& GXMatrix<T>::getRow(const unsigned& row)const{
    return this->mat[row];
}

template <typename T> std::vector<T>& GXMatrix<T>::getRow(const unsigned& row){
    return this->mat[row];
}


  // Access the row and column sizes  
template <typename T> unsigned GXMatrix<T>::get_rows() const{
	return this->rows;
} 
template <typename T> unsigned GXMatrix<T>::get_cols() const{
	return this->cols;
}

template <typename T> std::vector<std::vector<T> > GXMatrix<T>::getMat(){
	return this->mat;
}

template <typename T> std::string GXMatrix<T>::toString(){
	std::string s = "";
        int row = 0;
	for(std::vector<T> i : this->mat){
		s += "\n" + std::to_string(row) + ": ";
		for(T j : i) s += std::to_string(j) + "\t";
                row++;
	}
	return s;
}


#endif
