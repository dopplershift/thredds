package thredds.server.jupyter;

import javax.json.Json;
import javax.json.JsonObject;
import java.time.LocalDateTime;
import java.util.UUID;

class JupyterMessage {
    public enum MessageType {
        connect_request, connect_reply, status, stream,
        execute_request, execute_reply, execute_result, execute_input;
    };

    String session, user, msg_id, version;
    JsonObject metadata, content;
    LocalDateTime dateTime;
    MessageType type;
    JupyterMessage parent;

    public JupyterMessage(String session, String user, MessageType type) {
        this.user = user;
        this.type = type;
        this.session = session;

        this.metadata = null;
        this.content = null;
        this.parent = null;
        this.msg_id = UUID.randomUUID().toString();
        this.dateTime = LocalDateTime.now();
        this.version = "5.0";
    }

    public void addParent(JupyterMessage parent) {
        this.parent = parent;
    }

    public void addMetadata(JsonObject metadata) {
        this.metadata = metadata;
    }

    public void addContent(JsonObject content) {
        this.content = content;
    }

    public JsonObject toJson() {
       return Json.createObjectBuilder()
                .add("msg_id", msg_id)
                .add("username", user)
                .add("session", session)
                .add("msg_type", type.toString())
                .add("date", dateTime.toString())
                .add("version", version).build();
    }

    static public JupyterMessage fromJson(JsonObject json) {
        JupyterMessage ret = new JupyterMessage(json.getString("session"),
                json.getString("username"),
                MessageType.valueOf(json.getString("msg_type")));
        ret.msg_id = json.getString("msg_id");
        ret.dateTime = LocalDateTime.parse(json.getString("date"));
        ret.version = json.getString("version");
        return ret;
    }
}