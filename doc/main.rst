main module
===========


Call functions from PLinterfit library. An example can be as follows. 

    ::
    
    
       from PLinterafit.analysis import *
    
    
       # Press the green button in the gutter to run the script.
       if __name__ == '__main__':
           
           
           # Use this to fit 1 Voigt function to your data
           run_one_voigt_fitting()
           
           # Use this to fit 5 Voigt function to your data
           # run_five_voigt_fitting()
       
If you want to use your own data for the Voigt fitting, please change 'analysis' module according to your needs. 

.. automodule:: main
   :members:
   :undoc-members:
   :show-inheritance:
