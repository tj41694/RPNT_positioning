/* This is just for test! ---------------------------------------------------------------------------
* 2017/10/18	read O file
--------------------------------------------------------------------------------------------------- */
#include "BaseFunction/basefunction.h"
#include "GNSS/gnss_pro.h"

int main(int argc,char* argv[]) {
    int prd_Days=1;
    str2int(argv[2],prd_Days);

    all_option_c option;
    if (!option.readpostopt(argv[1])) return 0;
    for (int i=0; i<prd_Days; i++) {
        gnss_postsvr_c post;
        option.updatefile(i);

        if (!post.postsvrini(&option)) break;
       cout << "Day:" << setw(3) << i+1 << " start!\n";

        post.Post_Position_Epoch();
    }

    return 0;
}
