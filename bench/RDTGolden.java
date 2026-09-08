// A line-oriented adapter around the unmodified public RDT API.
import com.bioinceptionlabs.reactionblast.api.RDT;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.util.Base64;

public class RDTGolden {
    public static void main(String[] args) throws Exception {
        var cpu = (com.sun.management.OperatingSystemMXBean)
            ManagementFactory.getOperatingSystemMXBean();
        var input = new BufferedReader(new InputStreamReader(System.in, StandardCharsets.UTF_8));
        System.out.println("RDT_READY");
        String reaction;
        while ((reaction = input.readLine()) != null) {
            long wallStart = System.nanoTime(), cpuStart = cpu.getProcessCpuTime();
            String status, value;
            try {
                value = RDT.map(reaction).getMappedSmiles();
                status = "OK";
            } catch (Exception error) {
                value = error.toString();
                status = "ERROR";
            }
            double wall = (System.nanoTime() - wallStart) / 1e9;
            double usedCpu = (cpu.getProcessCpuTime() - cpuStart) / 1e9;
            String encoded = Base64.getEncoder().encodeToString(value.getBytes(StandardCharsets.UTF_8));
            System.out.println("RDT_RESULT\t" + status + "\t" + wall + "\t" + usedCpu + "\t" + encoded);
            System.out.flush();
        }
    }
}
