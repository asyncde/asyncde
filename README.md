### AsyncDE

AsyncDE is a high performance C++ implementation of the Asynchronous Differential Evolution algorithm, based on ideas from
<ul>
  <li>Zhabitskaya, E., Zhabitsky, M. (2013). Asynchronous Differential Evolution with Restart.</li>
 In: Dimov, I., Faragó, I., Vulkov, L. (eds) Numerical Analysis and Its Applications. NAA 2012. Lecture Notes in Computer Science, vol 8236. Springer,
 Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-41515-9_64
  <li>Zhabitsky, M., Zhabitskaya, E. 2013. Asynchronous differential evolution with adaptive correlation matrix.</li>
 In Proceedings of the 15th annual conference on Genetic and evolutionary computation (GECCO '13). Association for Computing Machinery, New York, NY, USA, 455–462.
 https://doi.org/10.1145/2463372.2463428
</ul>

The main algorithm &ndash; Asynchronous differential evolution with adaptive correlation matrix (ADE-ACM) &ndash; performs derivative-free optimization of multidimensional real-valued functions.
The algorithm has a simple structure and doesn't require tuning of its internal parameters by users. While performing the optimization, ADE-ACM self-adapts crossover rate and scale factor, and inflates the population size if needed.
The AsyncDE library also implements classic DE, JADE and jDE algorithms.

Feel free to contact authors with questions about AsyncDE by e-mail asyncde@gmail.com!

### Usage
AsyncDE requires C++11 compatible compiler (e.g. gcc or clang with development options) and CMake to compile. In addition, the library implements optional MPI and OMP variants for parallel calls of an optimized function within the Master/Workers scheme.

<!--
**asyncde/asyncde** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
