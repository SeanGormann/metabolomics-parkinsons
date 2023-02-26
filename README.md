# metabolomics-parkinsons
In this R project, I analyse data collected by Hu et al. (doi: 10.3389/fnmol.2020.00080). Specifically, I'm looking at the metabolomic profile of patients with parkinsons disease and age matched healthy controls. This data was obtained through LC-Mass Spec on plasma samples. In doing so, I aim to identify potential parkinson's biomarkers and gain information into the disease itself. 

---

After some standard data cleaning and preperation, sPLS-DA was preformed on the dataset:

![sPLSDA_metabolome](https://user-images.githubusercontent.com/100109163/221407495-260dfe5a-a9ca-4f5a-9a2e-ce451a4d561c.jpeg)


Through this, variables important in the projection (VIPs) were determined and used for further analysis. They can be visualised here:

![piechart_metabolomics](https://user-images.githubusercontent.com/100109163/221407568-47a1b1ea-56d0-4640-83c6-3952c8f569c8.jpeg)


And here: 

![wordcloud_metabolome](https://user-images.githubusercontent.com/100109163/221407598-33f76189-35ee-4d5e-9439-a7828108e9f0.jpeg)


---

The most significant VIPs were also identified. These were subjected to further analysis using the MetaboAnalyst software. Hypergeometric t-test were used to calculate pathway impacts. Files to be uploaded soon.

