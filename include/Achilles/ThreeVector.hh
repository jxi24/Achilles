#ifndef THREEVECTOR_HH
#define THREEVECTOR_HH

#include <array>
#include <cmath>
#include <iosfwd>

#include "spdlog/fmt/ostr.h"

namespace achilles {

/// @brief ThreeVector is a container for dealing with three vectors
///
/// The ThreeVector class provides an easy to use container to handle three
/// component vectors, such as position and three-momentums
class ThreeVector {
    using RotMat = std::array<double, 9>;

    public:
        /// @name Constructors and Destructors
        ///@{

        /// Create an empty ThreeVector object
        ThreeVector() noexcept : vec({{0, 0, 0}}) {}
        /// Create a ThreeVector object with values given by p
        ///@param p: A std::array<double, 3> containing the values for the vector
        ThreeVector(std::array<double, 3> p) noexcept : vec(p) {}
        /// Create a ThreeVector object with values pX, pY, and pZ
        ///@param pX: The x value of the vector
        ///@param pY: The y value of the vector
        ///@param pZ: The z value of the vector
        ThreeVector(double pX, double pY, double pZ) noexcept : vec({{pX, pY, pZ}}) {}
        /// Create a copy of a ThreeVector object
        ///@param other: The vector to be copied
        ThreeVector(const ThreeVector& other) noexcept = default;
        /// Move a ThreeVector object to another
        ///@param other: The vector to be moved
        ThreeVector(ThreeVector&& other) noexcept = default; 

        /// Assign a vector to another vector object
        ///@param other: three vector to be assigned
        ///@return ThreeVector: assigned new vector
        ThreeVector& operator=(const ThreeVector&) noexcept = default;
        ThreeVector& operator=(ThreeVector&&) noexcept = default;

        /// Default destructor
        ~ThreeVector() = default;
        ///@}

        /// @name Setters
        /// @{
        /// These functions provide access to setting the parameters of the ThreeVector

        /// Set position variable based on pased in array
        ///@param position: position to be stored
        void SetXYZ(const std::array<double, 3>& position) noexcept {
            vec = position;
        }

        /// Set position variable based on pased in array
        ///@param x: position in the x-direction to be stored
        ///@param y: position in the y-direction to be stored
        ///@param z: position in the z-direction to be stored
        void SetXYZ(const double& x, const double& y, const double& z) noexcept {
            vec = std::array<double, 3>{{x, y, z}};
        }

        /// Set position variable based on pased in array
        ///@param momentum: three momentum to be stored
        void SetPxPyPz(const std::array<double, 3> momentum) noexcept {
            vec = momentum;
        }

        /// Set momentum variable based on pased in array
        ///@param pX: momentum in the x-direction to be stored
        ///@param pY: momentum in the y-direction to be stored
        ///@param pZ: momentum in the z-direction to be stored
        void SetPxPyPz(const double& pX, const double& pY, const double& pZ) noexcept {
            vec = std::array<double, 3>{{pX, pY, pZ}};
        }

        /// Set only the position in the x-direction
        ///@param x: value to be stored
        void SetX(const double& x) noexcept {vec[0] = x;}

        /// Set only the position in the y-direction
        ///@param y: value to be stored
        void SetY(const double& y) noexcept {vec[1] = y;}

        /// Set only the position in the z-direction
        ///@param z: value to be stored
        void SetZ(const double& z) noexcept {vec[2] = z;}

        /// Set only the momentum in the x-direction
        ///@param pX: value to be stored
        void SetPx(const double& pX) noexcept {vec[0] = pX;}

        /// Set only the momentum in the y-direction
        ///@param pY: value to be stored
        void SetPy(const double& pY) noexcept {vec[1] = pY;}

        /// Set only the position in the z-direction
        ///@param pZ: value to be stored
        void SetPz(const double& pZ) noexcept {vec[2] = pZ;}
        ///@}

        /// @name Getters
        /// @{
        /// These functions provide get specific features from the ThreeVector object
        
        constexpr size_t Size() const noexcept { return 3; }

        /// Return the position as an array
        ///@return std::array<double, 3>: containing the position
        const std::array<double, 3>& Position() const noexcept { return vec; }

        /// Return the position in the x-direction
        ///@return double: containing the position in the x-direction
        const double& X() const noexcept { return vec[0]; }
        // double& X() noexcept { return vec[0]; }

        /// Return the position in the y-direction
        ///@return double: containing the position in the y-direction
        const double& Y() const noexcept { return vec[1]; }
        // double& Y() noexcept { return vec[1]; }

        /// Return the position in the z-direction
        ///@return double: containing the position in the z-direction
        const double& Z() const noexcept { return vec[2]; }
        // double& Z() noexcept { return vec[2]; }

        /// Return the momentum in the x-direction
        ///@return double: containing the momentum in the x-direction
        const double& Px() const noexcept { return vec[0]; }
        // double& Px() noexcept { return vec[0]; }

        /// Return the momentum in the y-direction
        ///@return double: containing the momentum in the y-direction
        const double& Py() const noexcept { return vec[1]; }
        // double& Py() noexcept { return vec[1]; }

        /// Return the momentum in the z-direction
        ///@return double: containing the momentum in the z-direction
        const double& Pz() const noexcept { return vec[2]; }
        // double& Pz() noexcept { return vec[2]; }

        /// Return the transverse momentum squared
        ///@return double: containing the transverse momentum squared
        double Pt2() const noexcept { return pow(vec[0], 2) + pow(vec[1], 2); }

        /// Return the transverse momentum
        ///@return double: containing the transverse momentum
        double Pt() const noexcept { return sqrt(Pt2()); }

        /// Return the three momentum squared
        ///@return double: containing the three momentum squared
        double P2() const noexcept { return (*this)*(*this); }

        /// Return the three momentum
        ///@return double: containing the three momentum
        double P() const noexcept { return sqrt(P2()); }

        /// Return the magnitude squared of the vector
        ///@return double: containing the magnitude squared of the vector
        double Magnitude2() const noexcept { return P2(); }

        /// Return the magnitude of the vector
        ///@return double: containing the magnitude of the vector
        double Magnitude() const noexcept { return P(); }

        /// Return the angle between the z-component and the transverse component
        ///@return double: containing the angle between the z and transverse components
        double Theta() const noexcept;

        /// Return the angle between the x and y components
        ///@return double: containing the angle between the x and y components
        double Phi() const noexcept;

        /// Rotate the three vector to the frame given by the 3 angles
        ///@param angles: The rotation angles
        ///@return FourVector: The vector in the corresponding frame
        ThreeVector Rotate(const std::array<double, 3>&) const noexcept;

        /// Rotate the three vector to the frame given by the 3 angles
        ///@param angles: The rotation matrix
        ///@return FourVector: The vector in the corresponding frame
        ThreeVector Rotate(const RotMat&) const noexcept;

        /// Obtain the rotation matrix to align the vector with a given axis
        ///@param axis: The axis to rotate to align with
        ///@return std::array<double, 9>: The rotation matrix to align the vector
        ///                               with the given axis
        RotMat Align(const ThreeVector&) const noexcept;

        /// Get the rotation matrix to align the vector with the z-axis
        ///@return std::array<double, 9>: The matrix needed to define the rotation
        RotMat AlignZ() const noexcept;
        ///@}

        /// @name Functions
        /// @{

        /// Calculate the dot product between two three vectors
        ///@param other: three vector to take dot product with
        ///@return double: the value of the dot product
        double Dot(const ThreeVector& other) const noexcept {
            return (*this) * other;
        }

        /// Calculate the cross product between two three vectors
        ///@param other: three vector to take cross product with
        ///@return ThreeVector: the resultant cross product
        ThreeVector Cross(const ThreeVector&) const noexcept;

        /// Return a unit vector in the direction of the vector
        ///@return ThreeVector: a unit vector in the direction of the vector
        ThreeVector Unit() const noexcept;

        /// Return a string representation of the vector
        ///@return std::string: a string representation of the vector
        std::string ToString() const noexcept;
        ///@}

        /// @name Operator Overloads
        /// @{
        /// Operator overloads of math functions for ease of use

        /// Add two three vectors together
        ///@param other: three vector to add to this one
        ///@return ThreeVector: The sum of the two three vectors
        ThreeVector& operator+=(const ThreeVector&) noexcept;

        /// Subtract two three vectors together
        ///@param other: three vector to subtract from this one
        ///@return ThreeVector: The difference of the two three vectors
        ThreeVector& operator-=(const ThreeVector&) noexcept;

        /// Scale a three vector by a constant
        ///@param scale: value to scale three vector by
        ///@return ThreeVector: The scaled three vector
        ThreeVector& operator*=(const double&) noexcept;

        /// Reduce the magnitude of a three vector by a constant
        ///@param scale: value to scale three vector by
        ///@return ThreeVector: The scaled three vector
        ThreeVector& operator/=(const double&);

        /// Calculate dot product of two vectors
        ///@param other: The other vector to take the dot product with
        ///@return double: The value of the dot product
        double operator*(const ThreeVector&) const noexcept;

        /// Negate a given vector
        ///@return ThreeVector: The negative of the input vector
        ThreeVector operator-() const noexcept;

        /// Unary plus operator
        ///@return ThreeVector: The input vector
        ThreeVector operator+() const noexcept;

        /// Scale a three vector by a constant
        ///@param scale: value to scale three vector by
        ///@return ThreeVector: The scaled three vector
        ThreeVector operator*(const double&) const noexcept;

        /// Reduce the magnitude of a three vector by a constant
        ///@param scale: value to scale three vector by
        ///@return ThreeVector: The scaled three vector
        ThreeVector operator/(const double&) const;

        /// Calculate difference of two vectors
        ///@param other: The other vector to take the difference with
        ///@return ThreeVector: The resulting difference vector
        ThreeVector operator-(const ThreeVector&) const noexcept;

        /// Calculate sum of two vectors
        ///@param other: The other vector to take the sum with
        ///@return ThreeVector: The resulting sum vector
        ThreeVector operator+(const ThreeVector&) const noexcept;

        // Comparison Operators

        /// Determine if two three vectors are equivalent
        ///@param other: The vector to compare against
        ///@return bool: True, if the vectors are equal otherwise False
        bool operator==(const ThreeVector&) const noexcept;

        /// Determine if two three vectors are not equivalent
        ///@param other: The vector to compare against
        ///@return bool: False, if the vectors are equal otherwise True
        bool operator!=(const ThreeVector& other) const noexcept {return !(*this == other);}

        // Access Operators

        /// Access a given index from the vector
        ///@param idx: Index to access
        ///@return double: The value of the vector at the given index
        double& operator[](const std::size_t& idx) {
            if(idx > 2) throw std::range_error("Max value is 2.");
            return vec[idx];
        }

        /// Access a given index from the vector
        ///@param idx: Index to access
        ///@return double: The value of the vector at the given index
        const double& operator[](const std::size_t& idx) const {
            if(idx > 2) throw std::range_error("Max value is 2.");
            return vec[idx];
        }
        /// @}

        ///@name Stream Operators
        /// @{
        /// Stream operators for writing to and reading from streams

        /// Write out a vector to an output stream
        ///@param ostr: Output stream to write to
        ///@param vec: The three vector to be written out
        template<typename OStream>
        friend OStream& operator<<(OStream& os, const ThreeVector& vec3) {
            os << "ThreeVector(" << vec3.Px() << ", " << vec3.Py() << ", " << vec3.Pz() << ")";
            return os;
        }

        /// Write in a vector to an input stream
        ///@param istr: Input stream to read from
        ///@param vec: The three vector to be read into
        friend std::istream& operator>>(std::istream&, ThreeVector&);
        /// @}

    private:
        std::array<double, 3> vec;
};

ThreeVector operator*(const double&, const ThreeVector&) noexcept;

}

#endif // end of include guard: THREEVECTOR_HH
