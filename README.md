# 6961Project
This is the project repo of ECE6961, we are currently working on CP-OFDM and ZP-OFDM <br>

Since there are 4 symbols, one symbol(real data) is 16 bits, and cp is 4 bits, which combines as 20 bits for 1 symbol with CP.<br>

For the channel model, we have h = 4x16 matrix, our method uses the function that does not consider the CP at all, which means we remove the CP during the channel <br>

Currently, we have done part 1 of the project, the CP-OFDM and ZP-OFDM. Just need to compare d_original vs d_result, you can see they are the same. Which means The CP and ZP are working correctly.<br>

Currently, we have done the project 1 both task 1 and task 2.<br>

We find the correct way to Trx, and calculate a_hat<br>

The Project_part2 is to run benchmark<br>
The project_part2_testbench is to run test_data<br>
