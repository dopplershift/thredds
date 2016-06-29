package thredds.server.jupyter;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.UUID;

import javax.json.Json;
import javax.json.JsonObject;

import org.zeromq.ZContext;
import org.zeromq.ZMQ;
import ucar.nc2.constants.CDM;

public class JupyterClient {

    private class JupyterSocket {
        final private String DELIM = "<IDS|MSG>";

        ZMQ.Socket socket;
        JupyterSocket(ZMQ.Socket socket, String URI) {
            this.socket = socket;
            socket.connect(URI);
        }

        void close() {
            socket.close();
        }

        private String readSocket() {
            return new String(socket.recv(), CDM.utf8Charset);
        }

        private JsonObject readSocketJson() {
            String resp = readSocket();
            if (resp.isEmpty())
                return null;
            return Json.createReader(new StringReader(resp)).readObject();
        }

        private boolean sendJson(JsonObject json) {
            if (json != null) {
                return socket.sendMore(json.toString());
            } else {
                return socket.sendMore("{}");
            }
        }

        private boolean sendJson(JupyterMessage msg) {
            if (msg != null) {
                return sendJson(msg.toJson());
            } else {
                return sendJson((JsonObject)null);
            }
        }

        JupyterMessage sendMessage(JupyterMessage msg) {
            socket.sendMore(DELIM);
            socket.sendMore("");
            sendJson(msg);
            sendJson(msg.parent);
            sendJson(msg.metadata);
            sendJson(msg.content);
            socket.send("");

            return receiveMessage();
        }

        JupyterMessage receiveMessage() {
            String delim = "";
            while (!delim.equals(DELIM)) {
                delim = readSocket();
            }

            String signature = readSocket();
            JupyterMessage resp = JupyterMessage.fromJson(readSocketJson());
            resp.addParent(JupyterMessage.fromJson(readSocketJson()));
            resp.addMetadata(readSocketJson());
            resp.addContent(readSocketJson());
            return resp;
        }
    }

    private class JupyterConnection {
        public String session;

        ZContext ctx;
        JupyterSocket cmdSocket, shellSocket, resultSocket;

        final int controlPort = 5000;

        public JupyterConnection() {
            ctx = new ZContext();
            cmdSocket = new JupyterSocket(ctx.createSocket(ZMQ.REQ),
                    "tcp://localhost:" + String.valueOf(controlPort));
            session = UUID.randomUUID().toString();
        }

        public void connect() {
            JupyterMessage control = new JupyterMessage(connection.session, "rmay",
                    JupyterMessage.MessageType.connect_request);
            JupyterMessage result = cmdSocket.sendMessage(control);
            shellSocket = new JupyterSocket(ctx.createSocket(ZMQ.DEALER),
                    "tcp://localhost:" + String.valueOf(result.content.getInt("shell")));

            // Must subscribe to a topic (even an empty one) to get messages
            ZMQ.Socket subSocket = ctx.createSocket(ZMQ.SUB);
            subSocket.subscribe("".getBytes());
            resultSocket = new JupyterSocket(subSocket,
                    "tcp://localhost:" + String.valueOf(result.content.getInt("iopub")));
        }

        public boolean execute(String code) {
            JupyterMessage exec = new JupyterMessage(connection.session, "rmay",
                    JupyterMessage.MessageType.execute_request);
            exec.content = Json.createObjectBuilder()
                    .add("code", code)
                    .add("silent", false)
                    .add("store_history", false)
                    .add("user_expressions", "")
                    .add("allow_stdin", false)
                    .add("stop_on_error", false)
                    .build();
            JupyterMessage result = shellSocket.sendMessage(exec);
            return result.content.getString("status").equals("ok");
        }

        public String getResult() {
            JupyterMessage msg;
            while ((msg = resultSocket.receiveMessage()) != null) {
                if (msg.type == JupyterMessage.MessageType.execute_result) {
                    return msg.content.getJsonObject("data").getString("text/plain");
                } else if (msg.type == JupyterMessage.MessageType.stream) {
                    return msg.content.getString("text");
                }
            }
            return "";
        }

        public void close() {
            if (resultSocket != null)
                resultSocket.close();
            if (shellSocket != null)
                shellSocket.close();
            cmdSocket.close();
            ctx.close();
        }
    }

    JupyterConnection connection;

    String connect(Path configPath) throws IOException {
        connection = new JupyterConnection();
        connection.connect();

        File codeFile = configPath.resolve("jupyter/test.py").toFile();
        FileInputStream fis = new FileInputStream(codeFile);
        byte[] codeBuff = new byte[(int) codeFile.length()];
        fis.read(codeBuff);
        if (connection.execute(new String(codeBuff, CDM.UTF8))) {
            String result = connection.getResult();
            return "Success! " + result;
        } else {
            return "Error executing code.";
        }
    }

    public void close() {
        connection.close();
    }
//        ProcessBuilder builder = new ProcessBuilder("/Users/rmay/miniconda3/envs/py35/bin/python",
//                "-m", "ipykernel", "--control", String.valueOf(controlPort));
//        builder.redirectErrorStream(true);
//        Process process = builder.start();
//
//        InputStream inp = process.getInputStream();
}
