<h1>Emulating Power Spectra of Dark Matter Halos using Neural Networks</h1>

<h2>Overveiw of the Project</h2>

<p>This project was part of a 'Senior Honours Project' course. The project aims to assertain the viability of using neural networks to model power spectra of dark matter halos using data taken from the Quijote Simulations Database.</p>

<h2>Abstract</h2>

<p>Three Tensorflow neural networks were optimised in the number of layers and nodes per layer,
training on data provided by the Quijote simulation. Attempting to emulate different regions of
the power spectra from five cosmological parameters of the ΛCDM model, in order to determine
which approach is most viable, if any. Three models were created: a low k model (ranging from
0.005 < k > 0.5), a high k model (ranging from 0.2 < k > 0.5) and a model also consisting of
a low k as an extra input for a high k output (low k ranging from 0.005 < k > 0.2, and high
k ranging from 0.2 < k > 0.5). The optimal number of layers were determined as 2, 2 and 3
with an optimal number of nodes of 600, 200, and 300 for the low k, the high k, and the low k
input, high k output models respectively. These models could each emulate the power spectra
accurately to a 5% error margin, with the Low k input, High k output model generating the most
accurate predictions when compared with the average values from the Quijote simulation’s P(k)
values. Further insights into whether the optimisation of the layers and nodes of each model were
the result of statistical randomness or not could be ascertained via more focused re-runs of the
optimisation functions with different seeds for training and testing samples, as well as a more
precise number for the number of nodes, in order to achieve less error in final model predictions.</p>

<h2>Results (Examples)</h2>

<p>Below are some example plots that were the results of the project (these can also be found in the <a href="/Senior Honours Project/Results/Dark_Matter_Halos_Neural_Networks.pdf">pdf</a>).</p>

<p>This first image shows the distribution of the data when split between training and testing data:</p>
<div align="center">
  <a href="/Senior Honours Project/Results/Images/all_models_cos_params.png"><img src="/Senior Honours Project/Results/Images/all_models_cos_params.jpg"></img></a>
</div>
</br>

<p>This image shows the difference between the true values and the predicted values for each neural network:</p>
<div align='center'>
  <a href='/Senior Honours Project/Results/Images/all_models_diffs.png'><img src='/Senior Honours Project/Results/Images/all_models_diffs.png'></img></a>
</div>
</br>

<p>This image shows a snippet of the difference between the true and predicted values for a specific model:</p>
<div align="center">
  <a href='/Senior Honours Project/Results/Images/lowk_highk_power_spectra.png'><img src="/Senior Honours Project/Results/Images/lowk_highk_power_spectra.png"></img></a>
</div>


<h2>Additional Info</h2>
<p>Due to issues with loss of access to the online portal containing my work (this is due to me completing my university course) there are no example notebook files showing full implementation of these functions, or files used for plotting the data.</p>
