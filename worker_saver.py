#
#   To simplify worker
#


class worker_container():
    def __init__(self, worker) -> None:

        self.StateName        = worker.StateName
        self.relative_error_tolerance = worker.relative_error_tolerance


        self.num_elements     = worker.num_elements
        self.num_labels       = worker.num_labels

        self.Noise            = worker.Noise
        if self.Noise == None:
            self.num_shots        = worker.num_shots

        self.label_list       = worker.label_list

        self.Nr               = worker.Nr
        self.meas_path        = worker.meas_path


        self.method           = worker.method
        
        self.Pj_method        = worker.Pj_method
        self.InitX            = worker.InitX

        self.measurement_list = worker.measurement_list

        self.InitErr          = worker.InitErr
        self.X0               = worker.X0
        self.Xk               = worker.Xk

        self.Err_relative_Xk  = worker.Err_relative_Xk
        self.Target_Err_Xk    = worker.Target_Err_Xk
        self.Target_Err_st    = worker.Target_Err_st
        self.fidelity_list    = worker.fidelity_list

        #self.target_relative_error_list = worker.target_relative_error_list      #  = []       

        self.step_Time        = worker.step_Time

        self.num_iterations   = worker.num_iterations
        self.iteration        = worker.iteration
        self.converged        = worker.converged

        self.method_dependent(worker)

    def method_dependent(self, worker):
        if self.method == 'RGD':
            self.Hermitian_List = worker.Hermitian_List
            self.EigV_pm_List   = worker.EigV_pm_List

            self.Ch_svd         = worker.Ch_svd

            self.u0             = worker.u0
            self.v0             = worker.v0
            self.s0             = worker.s0
            if self.InitX == 0:
                self.s0_choice  = worker.s0_choice

            self.uk             = worker.uk
            self.vk             = worker.vk
            self.ukh            = worker.ukh
            self.vkh            = worker.vkh

            self.sDiag          = worker.sDiag
            self.sDiag_list     = worker.sDiag_list
            self.Rec_Alpha      = worker.Rec_Alpha
            self.uGv_list       = worker.uGv_list

            self.coef           = worker.coef


        elif self.method == 'MiFGD':
            self.mu             = worker.mu

    def save_Run_Time(dt):
        self.RunT     =    dt           #   total Run time

    def Print_Info(self):
        print('     method = {}'.format(self.method))


if __name__ == "__main__":

    wc1 = worker_container(worker)        

