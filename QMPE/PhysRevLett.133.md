## Observing the Quantum Mpemba Effect in Quantum Simulations

Lata Kh. Joshi ,1,2,3 Johannes Franke ,1,4 Aniket Rath ,5 Filiberto Ares ,3 Sara Murciano ,6 Florian Kranzl ,1,4 [URL 🔗](https://orcid.org/0000-0001-6605-8027)

Rainer Blatt,1,4 Peter Zoller ,1,2 Benoît Vermersch ,1,2,5 Pasquale Calabrese,3,7 [URL 🔗](https://orcid.org/0000-0003-4014-1505)

Christian F. Roos ,1,4 and Manoj K. Joshi [URL 🔗](https://orcid.org/0000-0001-7121-8259)

1Institute for Quantum Optics and Quantum Information, Austrian Academy ofSciences, Technikerstraße 21a, 6020 Innsbruck, Austria

2University of Innsbruck, Institute for Theoretical Physics, Technikerstraße 21a, 6020 Innsbruck, Austria

3SISSA and INFN, via Bonomea 265, 34136 Trieste, Italy

4University of Innsbruck, Institute for Experimental Physics, Technikerstraße 25, 6020 Innsbruck, Austria

5Univiversit´e Grenoble Alpes, CNRS, LPMMC, 38000 Grenoble, France

6Walter Burke Institute for Theoretical Physics, and Department of Physics and IQIM, Caltech, Pasadena, California 91125, USA

7International Centre for Theoretical Physics (ICTP), Strada Costiera 11, 34151 Trieste, Italy

(Received 18 January 2024; accepted 29 March 2024; published 1 July 2024)

The nonequilibrium physics of many-body quantum systems harbors various unconventional phenom- ena. In this Letter, we experimentally investigate one of the most puzzling of these phenomena—the quantum Mpemba effect, where a tilted ferromagnet restores its symmetry more rapidly when it is farther from the symmetric state compared to when it is closer. We present the first experimental evidence of the occurrence of this effect in a trapped-ion quantum simulator. The symmetry breaking and restoration are monitored through entanglement asymmetry, probed via randomized measurements, and postprocessed using the classical shadows technique. Our findings are further substantiated by measuring the Frobenius distance between the experimental state and the stationary thermal symmetric theoretical state, offering direct evidence of subsystem thermalization.

[DOI: 10.1103/PhysRevLett.133.010402](https://doi.org/10.1103/PhysRevLett.133.010402)

Introduction.—When a system is brought out of equi- librium, it may exhibit phenomena that defy conventional wisdom. One particularly puzzling example is the Mpemba effect, first described as the phenomenon where hot water freezes faster than cold water [1], and then extended to a wide variety of systems [2–8]. An anomalous relaxation, reminiscent of the Mpemba effect, can manifest at zero temperature in isolated many-body quantum systems. Specifically, starting from a configuration that breaks a symmetry, its restoration can happen more rapidly when the initial state shows a greater degree of symmetry breaking [9]. This behavior, dubbed as quantum Mpemba effect (QMPE), is driven by entanglement and quantum fluctua- tions. The origin and ubiquity of the QMPE are active areas of research. For instance, in integrable systems, the con- ditions under which the QMPE occurs are now well understood [10]. However, the existence of this effect remains outstanding in generic quantum systems—such as nonintegrable and synthetic quantum many-body sys- tems. In this regard, we experimentally explore the QMPE in a chain of spins coupled via power-law decaying interactions. [URL 🔗](#page-0)

Present-day programmable quantum simulators, with their impeccable ability to create, manipulate, and analyze quantum states, provide us with excellent test beds to examine nonequilibrium dynamics in many-body quantum systems [11–14]. To address our objective of investigating [URL 🔗](#page-0)

the contribution ofeach subsystem, i.e., Q ¼ QA þQ¯A. The possible values q of the subsystem charge QA define the charge sectors. When a reduced density matrix is symmetric,

1,4,* [URL 🔗](#page-0)

the QMPE, we employ a trapped-ion quantum simulator

consisting of N ¼ 12 interacting spin-1=2 particles. The system is initialized into a product state where each spin points in the z direction and subsequently all spins are tilted by an angle θ from the z axis. The spin states designated by

θ ¼ 0 and π (we will refer to them as ferromagnetic states) are U(1) symmetric as they remain invariant under a rotation about the z axis; conversely for 0 < θ < π (tilted ferromagnets), the states explicitly break such symmetry. The tilted ferromagnetic state is then evolved with the engineered XY Hamiltonian [15–17] that also possesses the U(1) symmetry and determines the properties of time- evolved states. The U(1) symmetry implies the conserva- tion of the magnetization along the z axis, i.e., the conserved charge is defined as Q ¼ 1=2 P j σz j. In the thermodynamic limit, the reduced density matrix of a subsystem relaxes to a Gibbs ensemble [18–20], and the initially broken symmetry is restored [9]. The latter is a direct consequence of the Mermin-Wagner theorem which forbids spontaneous breaking of a continuous symmetry in one-dimensional systems at finite temperatures [21,22]. [URL 🔗](#page-0)

To set up an experimental indicator of symmetry break- ing, we consider a bipartition of the system as A ∪ ¯A. The chargeQthat generates the U(1) symmetry decomposes into


it is block diagonal in the eigenbasis of QA. Based on this property, the restoration of the symmetry, at the level of a subsystemA, is estimated by the “entanglement asymmetry” (EA), defined as [9], [URL 🔗](#page-0)

ΔSA ¼ log½Trðρ2 AÞ − log½Trðρ2 A;QÞ;

ð1Þ

where ρA;Q denotes the symmetrized counterpart of ρA, i.e.,

ρA;Q ¼ P q∈Z ΠqρAΠq. Each Πq is the projector onto the

charge sector q in the subsystemA. TheEA is a non-negative quantity that vanishes if and only if ρA is symmetric, that is

ρA ¼ ρA;Q [9]. The EA has proven to be an effective tool for investigating broken symmetries, not only in out-of-equi- librium many-body systems [10,23] but also in quantum field theories [24–26] and black hole physics [27].For present studies, both parts of the above equation can be estimated using the classical shadow formalism from ran- domized measurements [28,29]. [URL 🔗](#page-0)

In the current context, the QMPE manifests when an input state with a larger value of the EA, i.e., a greater degree of symmetry breaking, undergoes faster symmetry restoration than the one with a smaller value of the EA. Studying EA as a function of time, the QMPE is identified by the presence of a crossing between the EA curves for two states which are initialized at different degrees of symmetry breaking [9]. This route to symmetry restoration can also be mapped via a state distance, such as the Frobenius distance [30], which will be presented as a complementary approach to studying the QMPE. These studies performed on a quantum simulator not only allow us to investigate the QMPE in a nonintegrable system but also provide insights into the robustness of the effect under realistic physical conditions. In contrast to other quantum versions of the Mpemba effect, which require an external reservoir to drive the system out of equilibrium [31–37], our studies involve an isolated quantum system undergoing a unitary evolution of a pure quantum state. We note that during the preparation of this Letter, the Mpemba effect with a single ion coupled to an external reservoir has been reported [38,39]. These investigations complement our work while examining the QMPE in a different scenario. [URL 🔗](#page-0)

Our experiment marks the first observation of the QMPE in a many-body quantum system. We show that the pre- sence of integrability-breaking interactions and decohe- rence do not undermine its occurrence. This Letter presents experimental studies of symmetry restoration of a tilted ferromagnet for various scenarios. At first, we present our findings for nearly unitary dynamics of a tilted ferromag- net; the experimental results show that the QMPE occurs for the interacting XY spin chain. We also explore the QMPE under local disorders added to the interacting spin chain, and the scenario where tilted spins solely interact with the environmental noise. The former case reveals the QMPE for weak disorder strengths, however, when the disorder reaches a certain strength, it slows down the symmetry restoration, and subsequently weakens the

presence of the QMPE. In the latter case, we do not find the QMPE.

Setup.—We use a N ¼ 12 qubit trapped-ion quantum simulator to study the QMPE. A linear string of calcium ions is held in a Paul trap and laser-cooled to the motional ground state. The spin states are encoded into two long- lived electronic states, jS1=2;m ¼ 1=2i ≡ j1i ≡ j↓i and

jD5=2;m ¼

5=2i

≡

j0i ≡ j↑i and manipulated by a nar-

row-linewidth 729 nm laser. Further experimental details are given in Supplemental Material (SM) Secs. A, B, C, and D [40]. The experimental recipe for investigating the QMPE using randomized measurements is presented in Fig. 1(a). The spin state is initialized in a ferromagnetic [URL 🔗](#page-0)

product state j↓i⊗N and subsequently all spins are tilted by an angle θ using a laser beam that resonantly couples the two qubit states. The tilted state is then time evolved (quenched) with the desired Hamiltonian. After performing local random rotations U, site-resolved projective mea- surements are performed on the time-evolved state. The observables related to symmetry breaking and the state distance, i.e., the EA and the Frobenius distance for the subsystem of interest A, are estimated from the classical shadows. In Secs. E, F, and G ofSM [40] we provide details of the estimators. [URL 🔗](#page-0)

Symmetry restoration with a long-range spin-spin interaction.—The experimental data showing the QMPE are presented in Fig. 1(b). We choose three tilt angles: θ ¼ 0.2π, θ ¼ 0.33π, and θ ¼ 0.5π in increasing order of symmetry breaking at t ¼ 0, i.e., EA ¼ 0.64, 1.12, and 1.3 for NA ¼ 4, respectively. The tilt angle dependence of the EA at initial time t ¼ 0 is presented in the inset of Fig. 1(b) [URL 🔗](#page-0)

for subsystem size NA ¼ 4. The tilted ferromagnetic states are then quenched with a U(1) symmetric Hamiltonian engineered in our experiment. We engineer a power-law decaying XY interaction expressed as the Hamiltonian,

where σa i denotes the Pauli matrices for a ¼ x, y at lattice

site i ¼ 1;…;N. The realized interaction strength and range are J0 ≈ 560 rad=s and α ≈ 1, respectively.

Under the time evolution with the XY Hamiltonian, in the thermodynamic limit, the subsystem initialized into an asymmetric state is expected to attain the same symmetry of the Hamiltonian, i.e., the subsystem attains U(1) symmetry [9]. The symmetry restoration for subsystems of size [URL 🔗](#page-0)

NA ¼ 4 is monitored by measuring the EA at various times; see Fig. 1(b). The plotted EA is the average over all [URL 🔗](#page-0)

possible subsystems with size NA ¼ 4 from the central 6 lattice sites. We witness that the EA decreases and tends to zero for the three tilt angles, indicating that the subsystem state ρA attains U(1) symmetry. The striking feature is that the EA decays faster for a state at θ ¼ 0.5π than for states at [URL 🔗](#page-0)

θ ¼ 0.33π and 0.2π; i.e., in Fig. 1(b) the purple curve [URL 🔗](#page-0)


*FIG. 1. (a) Protocol to measure the EAwith a 1D array of 12 trapped ions using randomized measurements. After rotating the initially prepared ferromagnetic state by an angle θ, the system is evolved (quenched) with the Hamiltonian of interest for time t. Subsequently, randomized local rotations and projective measurements are performed. The obtained bitstrings (constructed from site-resolved fluorescence images) are then classically processed to estimate the EA. Under the quench dynamics with XY Hamiltonian (2): (b) experimentally measured EA, averaged over all possible subsystems of size NA ¼ 4 chosen from ion indices between 4 and 9, is plotted for angles θ ¼ 0.2π, 0.33π, and 0.5π. A fast restoration of symmetry is observed for the largest angle, confirming the quantum [URL 🔗](#page-0)*

*Mpemba effect (QMPE). The inset presents the EA at initial time t ¼ 0 as a function of the tilt angle θ. Solid curves denote theory including decoherences and symbols are experimental data. Error bars are calculated via the jackknife resampling method (see SM, Sec. H). (c) The measured EA is plotted for all connected 4-qubit subsystems of the whole chain (on the left) and for all possible 4-qubit subsystems (connected and disconnected) between ion numbers 4 to 9 (on the right). These plots show the independence of the QMPE on the subsystem. The x axis represents the subsystem index detailed in the main text.*

reaches zero before the green or the orange curve. Naively, one might assume that the state with the smallest initial EA would restore the symmetry at the earliest. The observed crossing in the EA curves confirms the QMPE. In Fig. 1(b), the solid lines are numerical simulations carried out for the experimental conditions, while including decoherence [URL 🔗](#page-0)

effects discussed in SM [40]. Notably, our system N ¼ 12 sufficiently captures the features of the QMPE. Decreasing the system size results in a finite size effect and affects the observation (see SM, Sec. I [40]). On the other hand, increasing the system does not substantiate the observation but increases the experimental complexity. [URL 🔗](#page-0)

We further examine the robustness of the QMPE with respect to the choice of subsystem A. For this, we evaluate

the EA for various subsystems with fixed size NA ¼ 4 and plot the experimental results in Fig. 1(c) for tilt angles [URL 🔗](#page-0)

θ ¼ 0.2π (green surface) and θ ¼ 0.5π (purple surface). In the left panel, we display the results for subsystems that are connected. We have considered all connected subsystems

of size NA ¼ 4 from one edge of the chain to the other edge.

There are nine of them which are ½1; 2; 3; 4; ½2; 3; 4; 5;

½3; 4; 5; 6;…, here the integers denote ion indices. Notably, the subsystems that lie at the edges of the ion string display the crossing of the EA curves at later times than those from the middle of the chain (bulk region). This is attributed to the effect of boundaries, as the subsystems from the edge regions are accompanied by a lesser number of neighboring particles than the ones in the bulk. This effect becomes more pronounced for smaller system sizes, which we further discuss in SM [40]. Furthermore, the choice of the subsystem is extensively studied by considering all connected and disconnected subsystems from the bulk. In the right panel of Fig. 1(c), we display the EA for the two tilt angles as a function of time and subsystem index. Here, [URL 🔗](#page-0)

we consider all subsystems of size NA ¼ 4 forming 15 subsystems out of the central 6 ions. In contrast to the left panel of Fig. 1(c), the crossing time of the EA for the two tilt angles does not depend upon the choice of the subsystem, which implies that the QMPE in the bulk region is robust against the choice of the subsystem. [URL 🔗](#page-0)


QMPE in the presence of disordered interactions.—A central theme in modern research is to understand how localization alters the relaxation and thermalization dynam- ics in spin systems [60]. The versatility of our experimental setup enables us to study this interesting problem through the lens of the EA and the QMPE. On the experimental side, we add transverse disorder terms to our XY Hamiltonian, thus engineering a disorder Hamiltonian of the type, [URL 🔗](#page-0)

The last term denotes static disorder terms realized in the experiment with site-dependent light shift laser beams [61]. The disorder fields are chosen randomly from a uniform [URL 🔗](#page-0)

distribution with disorder strengths hi ∈wð0;J0Þ. The Hamiltonian (3) is U(1) symmetric for each disorder realization and so can be exploited to study the QMPE. The tilted ferromagnetic states for two choices of tilt angles θ ¼ 0.2π and 0.5π are evolved under this Hamiltonian. We [URL 🔗](#page-0)

consider weak (w ¼ 6) and strong (w ¼ 14) disorder configurations. The experiment is repeated for 5 disorder sets of the weak and strong disorder cases (see SM [40] for the measured disorder values). In Fig. 2, we show the EA averaged over the 5 disorder realizations and all 15 subsystems of size 4 out of the central 6 ions. For the weak disorder, we observe the QMPE as the EA curves cross for the two tilt angles, indicating a faster drop of the [URL 🔗](#page-0)

EA for the input state with a larger EA at t ¼ 0. These results prove the robustness of the QMPE for weak disorders. Oppositely, for strong disorders, there is no crossing between the EA curves within the experimental

*FIG. 2. Time evolution of the EAwith the disorder Hamiltonian*

*(3). (a) In the presence of weak disorders, hi ∈6ð0;J0Þ,we observe crossing in the EA for the two tilt angles as a function of interaction time, thus revealing the QMPE in the measurements. [URL 🔗](#page-0)*

*(b) For the strong disorders, hi ∈14ð0;J0Þ, the QMPE is not visible in the measured data. The presented data are the averaged EAs for all 4 qubit subsystems out of the central 6 qubits. Solid curves denote theory and symbols are experimental data where the error bars are the standard deviation of the mean over the disorder realizations.*

time window. In this case, the strong disorders localize the interactions to single sites thus preventing the subsystem thermalization.

Relaxation of tilted ferromagnetic states under dephasing.—So far, we have studied the QMPE under near-unitary time evolution and investigated the QMPE

for tilted ferromagnets evolved under the XY and XY þ disorder Hamiltonian. Now we will explore a different scenario where the tilted ferromagnet evolves under a nonunitary evolution such as pure dephasing of the spins with the environmental noise. We turn the spin-spin interaction off and let the system freely evolve for time t. Because of the interaction with the fluctuating magnetic field in the laboratory and the phase noise of the laser beam (see SM [40]), dephasing takes place in the experiment. Under this scenario, the initially asymmetric state is expected to symmetrize. [URL 🔗](#page-0)

The experimental data are presented in Fig. 3 and the results are corroborated by numerical simulations. The measurement protocol is similar to the one presented in Fig. 1(a) except that the states are time evolved under the ambient experimental noise. In this case, we observe the restoration of symmetry, however, the crossing of the EA for various tilt angles is not observed. The experimental data and numerical simulations, see inset, imply that pure dephasing evolution does not show the QMPE while undergoing symmetry restoration of initially asymmetric states. In Fig. 3, solid lines correspond to numerical simulations, which are carried out for the experimentally measured dephasing rates listed in Sec. C of SM. [URL 🔗](#page-0)

Frobenius distance to observe QMPE.—In a comple- mentary approach, the route to symmetrization can also be probed by measuring the distance between two states—the

experimental time-evolved state ρAðtÞ and the diagonal

*FIG. 3. Symmetry restoration under a pure dephasing evolution*

*for three tilted ferromagnetic states: θ ¼ 0.2π, 0.33π, and 0.5π. Solid lines correspond to numerical simulations and data points are the experimental results averaged over 15 subsystems of*

*NA ¼ 4 generated from the central 6 ions and the error bars are calculated through jackknife resampling (see SM, Sec. H). Inset: simulation results are presented for late times. A reduction of the EA for all three tilt angles implies symmetry restoration, however, the QMPE is absent for the pure dephasing case since the EA curves do not cross.*


*FIG. 4. The Frobenius distance between the experimental time-evolved state and the theoretical equilibrium state is shown*

*for three tilt angles and subsystem size*

*NA ¼ 4. The state that is*

*initially farthest from its equilibrium state attains the latter at the earliest time; confirming the occurrence of the QMPE in quantum simulations via Frobenius distance measures. Solid lines are numerical simulations and the symbols are experi- mental data averaged over subsystems from the central region. The error bars are calculated via the jackknife resampling method (see SM, Sec. H).*

ensemble ρDE A [62]. The diagonal ensemble describes the average behavior of any observable at long times. It is a mixed state that is diagonal in the eigenbasis of the quenching XY Hamiltonian and, therefore, is U(1) sym- metric. According to the eigenstate thermalization hypoth- esis, for large systems, this ensemble is equivalent to a Gibbs ensemble [18–20]. From the data obtained in the randomized measurements, we evaluate the Frobenius [URL 🔗](#page-0)

distance [30] between ρAðtÞ and ρDE A (see more details in Sec. G of SM [40]). [URL 🔗](#page-0)

The results are presented in Fig. 4 for three tilt angles [URL 🔗](#page-0)

θ ¼ 0.2π, 0.33π, and 0.5π. From the present analysis, we note that the states which start farther away from the corresponding diagonal ensemble attain it earlier than those which start comparably closer. As we see, the state at an angle 0.5π begins at the largest distance and relaxes to the diagonal ensemble the fastest. Such observation is the manifestation of the QMPE via Frobenius distance mea- surements. Our findings constitute the first experimental evidence that the subsystem as a whole attains a stationary Gibbs ensemble. Although there exist previous experimen- tal studies on the thermalization of isolated quantum systems [63–68], they predominantly focused on specific local observables. [URL 🔗](#page-0)

Conclusions and outlook.—In this Letter, we have presented the first experimental demonstration of the QMPE. Notably, this phenomenon manifests itself at times much shorter than those relative to finite-size effects, such as revivals. This separation of time scales enables us to observe the QMPE distinctly in our experimental setup,

comprising a modest number of qubits (N ¼ 12), even in the presence of decoherence and disorder. Our investigation employs two distinct quantities, the entanglement asym- metry and the Frobenius distance, establishing that the

QMPE is not specific to a particular observable. This versatility paves the way for further investigation employ- ing a variety of theoretical and experimental tools.

Our focus in this study has been on a specific class of initial states—the tilted ferromagnets. However, our find- ings motivate future studies with other configurations to determine the experimental conditions for the occurrence of the QMPE in ergodic systems. While the role of entangle- ment in driving thermalization in isolated quantum systems is well understood, it will be important to elucidate the interplay between entanglement and quantum fluctuations for the observation of the QMPE in generic quantum systems. The measurement protocol used in this work is easily applicable to other platforms with individual control and readout capabilities such as arrays of atoms in optical lattices, Rydberg systems, or superconducting qubits [17,69–72]. This spotlights tantalizing opportunities for the experimental investigation of the QMPE, especially in higher dimensional systems where a richer phenomenology is expected due to the possibility of spontaneous symmetry breaking at finite temperatures. Our experimental findings may stimulate further theoretical investigations about the role played by dissipation and disorder in the dynamical restoration of symmetry. We believe that further studies on the QMPEwill provide new protocols for faster preparation of thermal states (albeit at the subsystem level). Such states can then be used as input states for quantum simulation experiments. The experimentally accessible quantities dis- cussed in the present manuscript, such as the EA and the Frobenius distance, will be valuable tools to assess the quality of the state preparation. [URL 🔗](#page-0)

The experimental team (J. F., F. K., M. K. J., C. R., R. B.) would like to acknowledge funding from the Institut für Quanteninformation GmbH and funding under Horizon Europe programme HORIZON-CL4-2022-QUANTUM-02- SGA via the project 101113690 (PASQuanS2.1). L. K. J., B. V., and P. Z. acknowledge funding from the Austrian Science Foundation [grant DOI: 10.55776/P35891] (FWF, P 32597 N). L. K. J. acknowledges financial support from the PNRR MUR Project No. PE0000023- NQST, and HPC and supercomputer facilities of the University of Innsbruck, where most of the numerical simulations were carried out. P. Z. acknowledges the funding from the Simons Collaboration on Ultra-Quantum Matter, which is a grant from the Simons Foundation (No. 651440). P. C. and F. A. are supported by ERC under Consolidator Grant No. 771536 (NEMO). S. M. acknowledges the support from the Caltech Institute for Quantum Information and Matter and the Walter Burke Institute for Theoretical Physics at Caltech. Work in Grenoble is funded by the French National Research Agency via the JCJC project QRand (No. ANR-20-CE47-0005), and via the research programs Plan France 2030 EPIQ (No. ANR-22-PETQ- 0007), QUBITAF (No. ANR-22-PETQ-0004), and HQI


(No. ANR-22-PNCQ-0002). A. R. acknowledges support by Laboratoire d’excellence LANEF in Grenoble (No. ANR-10-LABX-51-01) and from the Grenoble Nanoscience Foundation. In numerical simulations, we

have used the quantum toolbox

QuTiP [73]. [URL 🔗](#page-0)

*manoj.joshi@uibk.ac.at [URL 🔗](#page-0)

- [1] E. B. Mpemba and D. G. Osborne, Phys. Educ. 4, 312 (1969). [URL 🔗](https://doi.org/10.1088/0031-9120/4/3/312)

- [2] A. Lasanta, F. Vega Reyes, A. Prados, and A. Santos, Phys. Rev. Lett. 119, 148001 (2017). [URL 🔗](https://doi.org/10.1103/PhysRevLett.119.148001)

- [3] Z. Lu and O. Raz, (2017). Proc. Natl. Acad. Sci. U.S.A. 114, 5083 [URL 🔗](https://doi.org/10.1073/pnas.1701264114)

- [4] I. Klich, O. Raz, O. Hirschberg, and M. Vucelja, Phys. Rev. X 9, 021060 (2019). [URL 🔗](https://doi.org/10.1103/PhysRevX.9.021060)

- [5] A. Kumar and J. Bechhoefer, (2020). Nature (London) 584,64 [URL 🔗](https://doi.org/10.1038/s41586-020-2560-x)

- [6] J. Bechhoefer, A. Kumar, and R. Ch´etrite, 534 (2021). Nat. Rev. Phys. 3, [URL 🔗](https://doi.org/10.1038/s42254-021-00349-8)

- [7] A. Kumar, R. Ch´etrite, and J. Bechhoefer, Sci. U.S.A. 119, e2118484119 (2022). Proc. Natl. Acad. [URL 🔗](https://doi.org/10.1073/pnas.2118484119)

- [8] G. Teza, R. Yaacoby, and O. Raz, Phys. Rev. Lett. 131, 017101 (2023). [URL 🔗](https://doi.org/10.1103/PhysRevLett.131.017101)

- [9] F. Ares, S. Murciano, and P. Calabrese, Nat. Commun. 14, 2036 (2023). [URL 🔗](https://doi.org/10.1038/s41467-023-37747-8)

- [10] C. Rylands, K. Klobas, F. Ares, P. Calabrese, S. Murciano, and B. Bertini, arXiv:2310.04419. [URL 🔗](https://doi.org/10.1038/s41467-023-37747-8)

- [11] J. Zhang, G. Pagano, P. W. Hess, A. Kyprianidis, P. Becker, H. Kaplan, A. V. Gorshkov, Z.-X. Gong, and C. Monroe, Nature (London) 551, 601 (2017). [URL 🔗](https://arXiv.org/abs/2310.04419)

- [12] M. K. Joshi, F. Kranzl, A. Schuckert, I. Lovas, C. Maier, R. Blatt, M. Knap, and C. F. Roos, Science 376, 720 (2022). [URL 🔗](https://doi.org/10.1038/nature24654)

- [13] P. Scholl, M. Schuler, H. J. Williams, A. A. Eberharter, D. Barredo, K.-N. Schymik, V. Lienhard, L.-P. Henry, T. C. Lang, T. Lahaye, A. M. Läuchli, and A. Browaeys, Nature (London) 595, 233 (2021). [URL 🔗](https://doi.org/10.1126/science.abk2400)

- [14] J. W. Britton, B. C. Sawyer, A. C. Keith, C.-C. J. Wang, J. K. Freericks, H. Uys, M. J. Biercuk, and J. J. Bollinger, (London) 484, 489 (2012). Nature [URL 🔗](https://doi.org/10.1038/s41586-021-03585-1)

- [15] D. Porras and J. I. Cirac, Phys. Rev. Lett. 92, 207901 (2004). [URL 🔗](https://doi.org/10.1103/PhysRevLett.92.207901)

- [16] P. Jurcevic, B. P. Lanyon, P. Hauke, C. Hempel, P. Zoller, R. Blatt, and C. F. Roos, Nature (London) 511, 202 (2014). [URL 🔗](https://doi.org/10.1103/PhysRevLett.92.207901)

- [17] C. Monroe, W. C. Campbell, L.-M. Duan, Z.-X. Gong, A. V. Phys. 93, Gorshkov, P. W. Hess, R. Islam, K. Kim, N. M. Linke, G. Pagano, P. Richerme, C. Senko, and N. Y. Yao, Rev. Mod. 025001 (2021). [URL 🔗](https://doi.org/10.1038/nature13461)

- [18] J. M. Deutsch, Phys. Rev. A 43, 2046 (1991). [URL 🔗](https://doi.org/10.1103/RevModPhys.93.025001)

- [19] M. Srednicki, Phys. Rev. E 50, 888 (1994). [URL 🔗](https://doi.org/10.1103/PhysRevE.50.888)

- [20] M. Rigol, V. Dunjko, and M. Olshanii, Nature (London) 452, 854 (2008). [URL 🔗](https://doi.org/10.1038/nature06838)

- [21] P. C. Hohenberg, Phys. Rev. 158, 383 (1967). [URL 🔗](https://doi.org/10.1103/PhysRev.158.383)

- [22] N. D. Mermin and H. Wagner, Phys. Rev. Lett. 17, 1133 (1966). [URL 🔗](https://doi.org/10.1103/PhysRevLett.17.1133)

- [23] B. J. J. Khor, D. M. Kürkçüoglu, T. J. Hobbs, G. N. Perdue, and I. Klich, arXiv:2312.08601.

- [24] L. Capizzi and M. Mazzoni, J. High Energy Phys. 12 (2023) 144. [URL 🔗](https://doi.org/10.1007/JHEP12(2023)144)

- [25] M. Chen and H.-H. Chen, arXiv:2310.15480. [URL 🔗](https://arXiv.org/abs/2310.15480)

- [26] L. Capizzi and V. Vitale, arXiv:2310.01962. [URL 🔗](https://arXiv.org/abs/2310.01962)

- [27] F. Ares, S. Murciano, L. Piroli, and P. Calabrese, arXiv: 2311.12683. [URL 🔗](https://arXiv.org/abs/2311.12683)

- [28] A. Elben, S. T. Flammia, H.-Y. Huang, R. Kueng, J. Preskill, B. Vermersch, and P. Zoller, Nat. Rev. Phys. 5, 9 (2023). [URL 🔗](https://arXiv.org/abs/2311.12683)

- [29] H.-Y. Huang, R. Kueng, and J. Preskill, (2020). Nat. Phys. 16, 1050 [URL 🔗](https://doi.org/10.1038/s42254-022-00535-2)

- [30] M. Fagotti and F. H. L. Essler, Phys. Rev. B 87, 245107 (2013). [URL 🔗](https://doi.org/10.1103/PhysRevB.87.245107)

- [31] A. Nava and M. Fabrizio, Phys. Rev. B 100, 125102 (2019). [URL 🔗](https://doi.org/10.1103/PhysRevB.100.125102)

- [32] F. Carollo, A. Lasanta, and I. Lesanovsky, 127, 060401 (2021). Phys. Rev. Lett. [URL 🔗](https://doi.org/10.1103/PhysRevB.100.125102)

- [33] S. K. Manikandan, Phys. Rev. Res. 3, 043108 (2021). [URL 🔗](https://doi.org/10.1103/PhysRevResearch.3.043108)

- [34] S. Kochsiek, F. Carollo, and I. Lesanovsky, Phys. Rev. A 106, 012207 (2022). [URL 🔗](https://doi.org/10.1103/PhysRevA.106.012207)

- [35] F. Ivander, N. Anto-Sztrikacs, and D. Segal, Phys. Rev. E 108, 014130 (2023). [URL 🔗](https://doi.org/10.1103/PhysRevE.108.014130)

- [36] A. K. Chatterjee, S. Takada, and H. Hayakawa, Lett. 131, 080402 (2023). Phys. Rev. [URL 🔗](https://doi.org/10.1103/PhysRevE.108.014130)

- [37] A. K. Chatterjee, S. Takada, and H. Hayakawa, arXiv: 2311.01347. [URL 🔗](https://arXiv.org/abs/2311.01347)

- [38] S. A. Shapira, Y. Shapira, J. Markov, G. Teza, N. Akerman, O. Raz, and R. Ozeri, following Letter, Phys. Rev. Lett. 133, 010403 (2024). [URL 🔗](https://arXiv.org/abs/2311.01347)

- [39] J. Zhang, G. Xia, C.-W. Wu, T. Chen, Q. Zhang, Y. Xie, W.-B. Su, W. Wu, C.-W. Qiu, P.-x. Chen et al., arXiv: 2401.15951. [URL 🔗](https://doi.org/10.1103/PhysRevLett.133.010403)

- [40] See Supplemental Material at http://link.aps.org/ supplemental/10.1103/PhysRevLett.133.010402 for infor- mation regarding experimental details, error analysis, and data analysis, which includes Refs. [41–59]. [URL 🔗](http://link.aps.org/supplemental/10.1103/PhysRevLett.133.010402)

- [41] I. Arrazola, J. S. Pedernales, L. Lamata, and E. Solano, Sci. Rep. 6, 30534 (2016). [URL 🔗](https://doi.org/10.1038/srep30534)

- [42] P. Cieśliński, S. Imai, J. Dziewior, O. Gühne, L. Knips, W. Laskowski, J. Meinecke, T. Paterek, and T. V´ertesi, arXiv:2307.01251. [URL 🔗](https://doi.org/10.1038/srep30534)

- [43] J. Haah, A.W. Harrow, Z. Ji, X. Wu, and N. Yu, IEEE Trans. Inf. Theory 63, 5628 (2017). [URL 🔗](https://doi.org/10.1109/TIT.2017.2719044)

- [44] T. Brydges, A. Elben, P. Jurcevic, B. Vermersch, C. Maier, B. P. Lanyon, P. Zoller, R. Blatt, and C. F. Roos, Science 364, 260 (2019). [URL 🔗](https://doi.org/10.1109/TIT.2017.2719044)

- [45] A. Rath, R. van Bijnen, A. Elben, P. Zoller, and B. Vermersch, Phys. Rev. Lett. 127, 200503 (2021). [URL 🔗](https://doi.org/10.1126/science.aau4963)

- [46] K. J. Satzinger, Y.-J. Liu, A. Smith, C. Knapp et al., Science 374, 1237 (2021). [URL 🔗](https://doi.org/10.1126/science.abi8378)

- [47] Y. Zhou, P. Zeng, and Z. Liu, Phys. Rev. Lett. 125, 200502 (2020). [URL 🔗](https://doi.org/10.1103/PhysRevLett.125.200502)

- [48] A. Elben, R. Kueng, H.-Y. R. Huang, R. van Bijnen, C. Kokail, M. Dalmonte, P. Calabrese, B. Kraus, J. Preskill, P. Zoller, and B. Vermersch, (2020). Phys. Rev. Lett. 125, 200501

- [49] A. Neven, J. Carrasco, V. Vitale, C. Kokail, A. Elben, M. Dalmonte, P. Calabrese, P. Zoller, B. Vermersch, R. Kueng, and B. Kraus, npj Quantum Inf. 7, 152 (2021).

- [50] S. Imai, N. Wyderka, A. Ketterer, and O. Gühne, Lett. 126, 150501 (2021). Phys. Rev. [URL 🔗](https://doi.org/10.1038/s41534-021-00487-y)

- [51] C. Zhang, Y.-Y. Zhao, N. Wyderka, S. Imai, A. Ketterer, N.-N. Wang, K. Xu, K. Li, B.-H. Liu, Y.-F. Huang, C.-F. Li, G.-C. Guo, and O. Gühne, arXiv:2307.04382. [URL 🔗](https://doi.org/10.1103/PhysRevLett.126.150501)


- [52] A. Rath, C. Branciard, A. Minguzzi, and B. Vermersch, Phys. Rev. Lett. 127, 260501 (2021).

- [53] V. Vitale, A. Rath, P. Jurcevic, A. Elben, C. Branciard, and B. Vermersch, arXiv:2307.16882. [URL 🔗](https://doi.org/10.1103/PhysRevLett.127.260501)

- [54] L. K. Joshi, A. Elben, A. Vikram, B. Vermersch, V. Galitski, and P. Zoller, Phys. Rev. X 12, 011018 (2022). [URL 🔗](https://arXiv.org/abs/2307.16882)

- [55] F. Mezzadri, Not. Am. Math. Soc. 54, 592 (2007). [URL 🔗](https://doi.org/10.1103/PhysRevX.12.011018)

- [56] W. Hoeffding, Ann. Math. Stat. 19, 293 (1948). [URL 🔗](https://doi.org/10.1214/aoms/1177730196)

- [57] V. Vitale, A. Elben, R. Kueng, A. Neven, J. Carrasco, B. Kraus, P. Zoller, P. Calabrese, B. Vermersch, and M. Dalmonte, SciPost Phys. 12, 1 (2022). [URL 🔗](https://doi.org/10.1214/aoms/1177730196)

- [58] A. Rath, V. Vitale, S. Murciano, M. Votto, J. Dubail, R. Kueng, C. Branciard, P. Calabrese, and B. Vermersch, PRX Quantum 4, 010318 (2023). [URL 🔗](https://doi.org/10.21468/SciPostPhys.12.1.001)

- [59] C.-F. J. Wu, Ann. Stat. 14, 1261 (1986). [URL 🔗](https://doi.org/10.1214/aos/1176350142)

- [60] D. A. Abanin, E. Altman, I. Bloch, and M. Serbyn, Rev. Mod. Phys. 91, 021001 (2019). [URL 🔗](https://doi.org/10.1103/RevModPhys.91.021001)

- [61] C. Maier, T. Brydges, P. Jurcevic, N. Trautmann, C. Hempel, B. P. Lanyon, P. Hauke, R. Blatt, and C. F. Roos, Phys. Rev. Lett. 122, 050501 (2019). [URL 🔗](https://doi.org/10.1103/RevModPhys.91.021001)

- [62] M. Rigol, V. Dunjko, V. Yurovsky, and M. Olshanii, Phys. Rev. Lett. 98, 050405 (2007). [URL 🔗](https://doi.org/10.1103/PhysRevLett.98.050405)

- [63] M. Gring, M. Kuhnert, T. Langen, T. Kitagawa, B. Rauer, M. Schreitl, I. Mazets, D. A. Smith, E. Demler, and J. Schmiedmayer, Science 337, 1318 (2012). [URL 🔗](https://doi.org/10.1103/PhysRevLett.98.050405)

- [64] S. Trotzky, Y.-A. Chen, A. Flesch, I. P. McCulloch, U. Schollwöck, J. Eisert, and I. Bloch, Nat. Phys. 8, 325 (2012).

- [65] A. M. Kaufman, M. E. Tai, A. Lukin, M. Rispoli, R. Schittko, P. M. Preiss, and M. Greiner, Science 353, 794 (2016).

- [66] C. Neill et al., Nat. Phys. 12, 1037 (2016). [URL 🔗](https://doi.org/10.1038/nphys3830)

- [67] M. Ueda, Nat. Rev. Phys. 2, 669 (2020). [URL 🔗](https://doi.org/10.1038/s42254-020-0237-x)

- [68] Z.-Y. Zhou, G.-X. Su, J. C. Halimeh, R. Ott, H. Sun, P. Hauke, B. Yang, Z.-S. Yuan, J. Berges, and J.-W. Pan, Science 377, 311 (2022). [URL 🔗](https://doi.org/10.1038/s42254-020-0237-x)

- [69] A. J. Daley, I. Bloch, C. Kokail, S. Flannigan, N. Pearson, M. Troyer, and P. Zoller, Nature (London) 607, 667 (2022). [URL 🔗](https://doi.org/10.1126/science.abl6277)

- [70] I. Bloch, J. Dalibard, and W. Zwerger, Rev. Mod. Phys. 80, 885 (2008). [URL 🔗](https://doi.org/10.1103/RevModPhys.80.885)

- [71] S. J. Evered, D. Bluvstein, M. Kalinowski, S. Ebadi, T. Manovitz, H. Zhou, S. H. Li, A. A. Geim, T. T. Wang, N. Maskara, H. Levine, G. Semeghini, M. Greiner, V. Vuletić, and M. D. Lukin, Nature (London) 622, 268 (2023). [URL 🔗](https://doi.org/10.1103/RevModPhys.80.885)

- [72] J. C. Hoke, M. Ippoliti, E. Rosenberg et al., Nature (London) 622, 481 (2023). [URL 🔗](https://doi.org/10.1038/s41586-023-06505-7)

- [73] J. Johansson, P. Nation, and F. Nori, Comput. Phys. Commun. 184, 1234 (2013). [URL 🔗](https://doi.org/10.1016/j.cpc.2012.11.019)
